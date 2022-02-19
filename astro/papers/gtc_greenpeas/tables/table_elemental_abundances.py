import numpy as np
import src.specsiser as sr
import pylatex
from delete.data_printing import PdfPrinter, latex_labels, numberStringFormat
from pathlib import Path

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
refernList = obsData['file_information']['refName_list']
resultsFolder = Path(obsData['file_information']['results_folder'])
tables_folder = Path(obsData['file_information']['tables_folder'])

ext = 'BR'
cycle = 'it3'

sample_dict = {}
for i, obj in enumerate(objList[:3]):

    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
    outputPiDb = objFolder / f'{obj}_{ext}_pI_fitting_{cycle}.db'
    outputDirectHIIchemDb = objFolder / f'{obj}_{ext}_Direct-Teff-logU_{cycle}.db'

    results_dict = sr.loadConfData(results_file, group_variables=False)
    section_label = f'Elemental_abundances_{cycle}'
    section_dict = results_dict[section_label]
    sample_dict[obj] = section_dict

# Generate the table
table_headers = ['Parameter'] + list(refernList[:3])

pdf = PdfPrinter()
pdf.create_pdfDoc(pdf_type=None)
pdf.pdf_insert_table(table_headers, table_format='c'*len(table_headers), addfinalLine=True)

param_list = list(sample_dict[objList[0]].keys())

latex_labels['O'] = r'$\nicefrac{O}{H}$'
latex_labels['N'] = r'$\nicefrac{N}{H}$'
latex_labels['S'] = r'$\nicefrac{S}{H}$'
latex_labels['NO'] = r'$\nicefrac{N}{O}$'
latex_labels['He'] = r'$\nicefrac{He}{H}$'

for j, param in enumerate(['O', 'N', 'S', 'He']):

    row_data = ['-'] * len(table_headers)
    row_data[0] = latex_labels[param]

    for i, obj in enumerate(objList[:3]):

        # Spherical fitting HII-CHI-mistry
        param_value = sample_dict[obj][param]
        if param == 'Teff':
            param_entry = r'${:.0f}\pm{:.0f}$'.format(param_value[0], param_value[1])
        else:
            param_entry = r'${:.2f}\pm{:.2f}$'.format(param_value[0], param_value[1])
        row_data[1 + i] = param_entry

    print(row_data)
    pdf.addTableRow(row_data)

pdf.table.add_hline()

for j, param in enumerate(['NO', 'ICF_SIV']):

    row_data = ['-'] * len(table_headers)
    row_data[0] = latex_labels[param]

    for i, obj in enumerate(objList[:3]):

        # Spherical fitting HII-CHI-mistry
        param_value = sample_dict[obj][param]
        if param == 'Teff':
            param_entry = r'${:.0f}\pm{:.0f}$'.format(param_value[0], param_value[1])
        else:
            param_entry = r'${:.2f}\pm{:.2f}$'.format(param_value[0], param_value[1])
        row_data[1 + i] = param_entry

    pdf.addTableRow(row_data)

pdf.table.add_hline()

table_file = tables_folder / f'elemental_abundances'
pdf.generate_pdf(table_file)

