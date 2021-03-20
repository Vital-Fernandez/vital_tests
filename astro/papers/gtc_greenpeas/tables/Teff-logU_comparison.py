import numpy as np
import src.specsiser as sr
import pylatex
from src.specsiser.data_printing import PdfPrinter, latex_labels, numberStringFormat
from pathlib import Path

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
refernList = obsData['file_information']['refName_list']
resultsFolder = Path(obsData['file_information']['results_folder'])
tables_folder = Path(obsData['file_information']['tables_folder'])


ext = 'BR'
cycle = 'it3'

table_headers = ['Galaxy',
                 'Parameter',
                 r'\makecell{\textsc{HII-CHI-mistry: } \\ Spherical}',
                 r'\makecell{\textsc{HII-CHI-mistry: } \\ Plane-Parallel}',
                 r'\makecell{Bayesian \textsc{HII-CHI-mistry: } \\ Plane-Parallel}',
                 r'\makecell{Bayesian \\ Direct Method + \textsc{HII-CHI-mistry} \\ Plane-Parallel}'
                 ]

pdf = PdfPrinter()
pdf.create_pdfDoc(pdf_type=None)
pdf.pdf_insert_table(table_headers, addfinalLine=True)

for i, obj in enumerate(objList[:3]):

    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
    outputPiDb = objFolder / f'{obj}_{ext}_pI_fitting_{cycle}.db'
    outputDirectHIIchemDb = objFolder / f'{obj}_{ext}_Direct-Teff-logU_{cycle}.db'

    results_dict = sr.loadConfData(results_file, group_variables=False)

    for j, param in enumerate(('Teff', 'logU')):

        row_data = ['-'] * len(table_headers)

        if j == 0:
            row_data[0] = pylatex.MultiRow(2, data=refernList[i])
        else:
            row_data[0] = ''

        row_data[1] = latex_labels[param]

        # Spherical fitting HII-CHI-mistry
        grid_file = 'C17_bb_Teff_30-90_sph.dat'
        section_label = f'HII_Tef_fit_{grid_file}_{cycle}_sed'
        param_value = results_dict[section_label][param]
        if param == 'Teff':
            param_entry = r'${:.0f}\pm{:.0f}$'.format(param_value[0], param_value[1])
        else:
            param_entry = r'${:.2f}\pm{:.2f}$'.format(param_value[0], param_value[1])
        row_data[2] = param_entry

        # Planar fitting HII-CHI-mistry
        grid_file = 'C17_bb_Teff_30-90_pp.dat'
        section_label = f'HII_Tef_fit_{grid_file}_{cycle}_sed'
        param_value = results_dict[section_label][param]
        if param == 'Teff':
            param_entry = r'${:.0f}\pm{:.0f}$'.format(param_value[0], param_value[1])
        else:
            param_entry = r'${:.2f}\pm{:.2f}$'.format(param_value[0], param_value[1])
        row_data[3] = param_entry

        # Bayesian HII-CHI-mistry
        model_fit = sr.load_MC_fitting(outputPiDb)
        trace = model_fit['trace'][param]
        param_value = (np.mean(trace), np.std(trace))
        if param == 'Teff':
            param_entry = r'${:.0f}\pm{:.0f}$'.format(param_value[0], param_value[1])
        else:
            param_entry = r'${:.2f}\pm{:.2f}$'.format(param_value[0], param_value[1])
        row_data[4] = param_entry

        # Bayesian Direct method + HII-CHI-mistry
        if outputDirectHIIchemDb.is_file():
            model_fit = sr.load_MC_fitting(outputDirectHIIchemDb)
            trace = model_fit['trace'][param]
            param_value = (np.mean(trace), np.std(trace))
            if param == 'Teff':
                param_entry = r'${:.0f}\pm{:.0f}$'.format(param_value[0], param_value[1])
            else:
                param_entry = r'${:.2f}\pm{:.2f}$'.format(param_value[0], param_value[1])
            row_data[5] = param_entry

        pdf.addTableRow(row_data)

pdf.table.add_hline()

table_file = tables_folder / f'photo-ionization_param_comparison'
pdf.generate_pdf(table_file)

