import src.specsiser as sr
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from src.specsiser.data_printing import PdfPrinter
import pyneb as pn

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
refName_list = obsData['file_information']['refName_list']
refName_list = ['GP030321', 'GP101157', 'GP121903', 'GP004054', 'GP113303', 'GP232539']


dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
tables_folder = Path(obsData['file_information']['tables_folder'])
fileList = obsData['file_information']['files_list']
idx_band = int(obsData['file_information']['band_flux'])
z_objs = obsData['sample_data']['z_array']
flux_norm = obsData['sample_data']['norm_flux']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

ext = 'BR'
cycle = 'it3'

table_labels = dict(CHI2_RED=r'$\chi^2$',
                    BST_MMET=r'$\left\langle z_{M}\right\rangle $',
                    BST_LMET=r'$\left\langle z_{L}\right\rangle $',
                    BSTLLAGE=r'$\left\langle log(t_{L}) \right\rangle $',
                    BSTLMAGE=r'$\left\langle log(t_{M}) \right\rangle $',
                    LOGMEBST=r'$log(M_{e})$',
                    LOGMCBST=r'$log(M_{c})$',
                    LOGPEBST=r'$log(M^{> 1Myr}_{e})$',
                    LOGPCBST=r'$log(M^{> 1Myr}_{c})$',
                    GEXTINCT=r'$c(H\beta)_{stellar}$',
                    GNEBULAR=r'$c(H\beta)_{nebular}$')
nebular_percetages = ['10.08', '27.57',	'23.28', '17.19', '10.92', '7.40']

# table_headers = ['Galaxy',
#                  r'$\chi^2$',
#                  r'$\left\langle z_{M}\right\rangle $-$\left\langle z_{L}\right\rangle $',
#                  r'$\left\langle log(t_{L}) \right\rangle $-$\left\langle log(t_{M}) \right\rangle $',
#                  r'$log(M_{e}) - log(M^{> 1Myr}_{e})$',
#                  r'$log(M_{c}) - log(M^{> 1Myr}_{c})$',
#                  r'$A_{V,\,stellar} - A_{V,\,nebular}$']

table_headers = ['Galaxy',
                 r'$\chi^2$',
                 r'$\left\langle z_{M}\right\rangle$',
                 r'$\left\langle log(t_{L}) \right\rangle$',
                 r'$log(M_{e})$',
                 r'$log(M_{c})$',
                 r'$A_{V,\,stellar}$']

table_headers2 = ['Galaxy',
                 r'$L_{nebular}\,(\%)$',
                 r'$\left\langle z_{L}\right\rangle $',
                 r'$\left\langle log(t_{M}) \right\rangle $',
                 r'$log(M^{> 1Myr}_{e})$',
                 r'$log(M^{> 1Myr}_{c})$',
                 r'$A_{V,\,nebular}$']

pdf = PdfPrinter()
pdf.create_pdfDoc(pdf_type=None)
pdf.pdf_insert_table(table_headers, addfinalLine=True)

for i, obj in enumerate(objList):

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'

    # Load the data
    results_dict = sr.loadConfData(results_file, group_variables=False)
    section_label = f'Large_FADO_fit'
    section_results = results_dict[section_label]

    row_data = [''] * len(table_headers)
    row_data[0] = refName_list[i]

    # Chi
    j = 1
    param_code = 'CHI2_RED'
    param_errCode = param_code.replace('BST', 'DEV')
    param_value = section_results[param_code]
    param_err = section_results[param_errCode] if param_errCode in section_results else None
    param_entry = r'${}$'.format(param_value)
    row_data[j] = param_entry

    # Metallicity
    j += 1
    param_code = ['BST_MMET', 'BST_LMET']
    param_entry = r'${:.3f}$'.format(section_results[param_code[0]])
    row_data[j] = param_entry

    # Age
    j += 1
    param_code = ['BSTLLAGE', 'BSTLMAGE']
    param_entry = r'${:.3f}$'.format(section_results[param_code[0]])
    row_data[j] = param_entry

    # Mass ever formed
    j += 1
    param_code = ['LOGMEBST', 'LOGPEBST']
    section_results[param_code[1]] = 0 if section_results[param_code[1]] < 0.0 else section_results[param_code[1]]
    param_entry = r'${:.3f}$'.format(section_results[param_code[0]])
    row_data[j] = param_entry

    # Mass current ever
    j += 1
    param_code = ['LOGMCBST', 'LOGPCBST']
    section_results[param_code[1]] = 0 if section_results[param_code[1]] < 0.0 else section_results[param_code[1]]
    param_entry = r'${:.3f}$'.format(section_results[param_code[0]])
    row_data[j] = param_entry

    # Mass current ever
    j += 1
    param_code = ['GEXTINCT', 'GNEBULAR']
    param_entry = r'${:.3f}$'.format(section_results[param_code[0]])
    row_data[j] = param_entry

    # Nebular fraction
    # j += 1
    # row_data[j] = r'${}\%$'.format(nebular_percetages[i])

    pdf.addTableRow(row_data, last_row=False)

pdf.table.add_hline()
pdf.addTableRow(table_headers2, last_row=False)
pdf.table.add_hline()

for i, obj in enumerate(objList):

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'

    # Load the data
    results_dict = sr.loadConfData(results_file, group_variables=False)
    section_label = f'Large_FADO_fit'
    section_results = results_dict[section_label]

    row_data = [''] * len(table_headers)
    row_data[0] = refName_list[i]

    # Chi
    j = 1
    # Nebular fraction
    # j += 1
    row_data[j] = r'${}\%$'.format(nebular_percetages[i])

    # Metallicity
    j += 1
    param_code = ['BST_MMET', 'BST_LMET']
    param_entry = r'${:.3f}$'.format(section_results[param_code[1]])
    row_data[j] = param_entry

    # Age
    j += 1
    param_code = ['BSTLLAGE', 'BSTLMAGE']
    param_entry = r'${:.3f}$'.format(section_results[param_code[1]])
    row_data[j] = param_entry

    # Mass ever formed
    j += 1
    param_code = ['LOGMEBST', 'LOGPEBST']
    section_results[param_code[1]] = 0 if section_results[param_code[1]] < 0.0 else section_results[param_code[1]]
    param_entry = r'${:.3f}$'.format(section_results[param_code[1]])
    row_data[j] = param_entry

    # Mass current ever
    j += 1
    param_code = ['LOGMCBST', 'LOGPCBST']
    section_results[param_code[1]] = 0 if section_results[param_code[1]] < 0.0 else section_results[param_code[1]]
    param_entry = r'${:.3f}$'.format(section_results[param_code[1]])
    row_data[j] = param_entry

    # Mass current ever
    j += 1
    param_code = ['GEXTINCT', 'GNEBULAR']
    param_entry = r'${:.3f}$'.format(section_results[param_code[1]])
    row_data[j] = param_entry

    pdf.addTableRow(row_data, last_row=False)

pdf.table.add_hline()

pdf.generate_pdf(tables_folder/'fado_fitting_results')



# for i, obj in enumerate(objList):
#
#     # Declare input files
#     objFolder = resultsFolder / f'{obj}'
#     results_file = objFolder / f'{obj}_{ext}_measurements.txt'
#
#     # Load the data
#     results_dict = sr.loadConfData(results_file, group_variables=False)
#     section_label = f'Large_FADO_fit'
#     section_results = results_dict[section_label]
#
#     row_data = [''] * len(table_headers)
#     row_data[0] = refName_list[i]
#
#     # Chi
#     j = 1
#     param_code = 'CHI2_RED'
#     param_errCode = param_code.replace('BST', 'DEV')
#     param_value = section_results[param_code]
#     param_err = section_results[param_errCode] if param_errCode in section_results else None
#     param_entry = r'${}$'.format(param_value)
#     row_data[j] = param_entry
#
#     # Metallicity
#     j += 1
#     param_code = ['BSTLLAGE', 'BST_LMET']
#     param_entry = r'${:.3f};\,{:.3f}$'.format(section_results[param_code[0]], section_results[param_code[1]])
#     row_data[j] = param_entry
#
#     # Age
#     j += 1
#     param_code = ['BSTLLAGE', 'BSTLMAGE']
#     param_entry = r'${:.3f};\,{:.3f}$'.format(section_results[param_code[0]], section_results[param_code[1]])
#     row_data[j] = param_entry
#
#     # Mass ever formed
#     j += 1
#     param_code = ['LOGMEBST', 'LOGPEBST']
#     section_results[param_code[1]] = 0 if section_results[param_code[1]] < 0.0 else section_results[param_code[1]]
#     param_entry = r'${:.3f};\,{:.3f}$'.format(section_results[param_code[0]], section_results[param_code[1]])
#     row_data[j] = param_entry
#
#     # Mass current ever
#     j += 1
#     param_code = ['LOGMCBST', 'LOGPCBST']
#     section_results[param_code[1]] = 0 if section_results[param_code[1]] < 0.0 else section_results[param_code[1]]
#     param_entry = r'${:.3f};\,{:.3f}$'.format(section_results[param_code[0]], section_results[param_code[1]])
#     row_data[j] = param_entry
#
#     # Mass current ever
#     j += 1
#     param_code = ['GEXTINCT', 'GNEBULAR']
#     param_entry = r'${:.3f};\,{:.3f}$'.format(section_results[param_code[0]], section_results[param_code[1]])
#     row_data[j] = param_entry
#
#     # Nebular fraction
#     j += 1
#     row_data[j] = r'${}\%$'.format(nebular_percetages[i])
#
#     pdf.addTableRow(row_data, last_row=False)
#
# pdf.table.add_hline()
# pdf.generate_pdf(tables_folder/'fado_fitting_results')


