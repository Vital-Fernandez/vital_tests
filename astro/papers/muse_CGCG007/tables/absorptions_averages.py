import numpy as np
from pathlib import Path

import lime
from lime.plots import PdfMaker
from lime import load_cfg

# Declare data and files location
obsData = load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
tablesFolder = Path(obsData['data_location']['tables_folder'])

abs_data = obsData['CGCG007_absorptions']
count_data = obsData['CGCG007_abs_count']

# Table headers
table_headers = ['Line', r'Relative absorption', 'Voxels']

# Pdf object
pdf = PdfMaker()
pdf.create_pdfDoc(pdf_type=None)
pdf.pdf_insert_table(table_headers, addfinalLine=True)

for param, value in abs_data.items():

    line = param.split('_abs')[0]
    ion, wave, latexLabel = lime.label_decomposition(line, scalar_output=True)
    absF, errF = value
    n_voxels = count_data[f'{line}_count']

    row = ['-'] * len(table_headers)

    row[0] = latexLabel
    row[1] = r'${} \pm {}$'.format(f'{absF:.2f}', f'{errF:.2f}') if line != 'H1_4861A' else '1.0'
    row[2] = str(int(n_voxels))
    print(row)

    pdf.addTableRow(row, last_row=False)

pdf.table.add_hline()
pdf.generate_pdf(tablesFolder/'absorptions')