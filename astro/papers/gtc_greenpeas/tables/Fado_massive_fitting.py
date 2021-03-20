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

table_headers = ['Galaxy'] + list(table_labels.values()) + ['Nebular contribution']


pdf = PdfPrinter()
pdf.create_pdfDoc(pdf_type='table')
pdf.pdf_insert_table(table_headers, addfinalLine=True)

for i, obj in enumerate(objList):

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'

    # Load the data
    results_dict = sr.loadConfData(results_file, group_variables=False)
    section_label = f'Large_FADO_fit'
    section_results = results_dict[section_label]

    print(refName_list[i])
    row_data = [refName_list[i]] + [''] * (len(table_labels) + 1)
    for j, item_pair in enumerate(table_labels.items()):

        param_code, param_label = item_pair
        param_errCode = param_code.replace('BST', 'DEV')
        param_value = section_results[param_code]
        param_err = section_results[param_errCode] if param_errCode in section_results else None

        # if param_code is None:
        #     param_entry = param_value
        # else:
        #     param_entry = r'${}\pm{}$'.format(param_value, param_err)

        if param_value == -999.0:
            param_entry = '-'
        else:
            if param_code in ['GEXTINCT', 'GNEBULAR']:
                rc = pn.RedCorr(R_V=RV, E_BV=param_value / RV, law=red_law)
                param_err = section_results['GEXTBDEV'] if param_code is 'GEXTINCT' else section_results['GNEBBDEV']
                cHbeta_value, cHbeta_err = rc.cHbetaFromEbv(param_value / RV), rc.cHbetaFromEbv(param_err / RV)
                param_entry = r'${:.2f}\pm{:.3f}$'.format(cHbeta_value, cHbeta_err)

            else:
                param_entry = r'${}$'.format(param_value)

        row_data[j+1] = param_entry

    # Add the nebular contribution
    row_data[-1] = r'${}\%$'.format(nebular_percetages[i])

    pdf.addTableRow(row_data, last_row=False)

pdf.table.add_hline()
pdf.generate_pdf(tables_folder/'fado_fitting_results')
