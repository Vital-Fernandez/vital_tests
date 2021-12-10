# import src.specsiser as sr
# from pathlib import Path
# import matplotlib.pyplot as plt
# import numpy as np
# from src.specsiser.data_printing import PdfPrinter
# import pyneb as pn

import numpy as np
from pathlib import Path
from lime.io import PdfMaker
from lime import load_cfg

# Declare data and files location
obsData = load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
tablesFolder = Path(obsData['data_location']['tables_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']


# Table headers
table_headers = ['Parameter', 'Value']

# Pdf object
pdf = PdfMaker()
pdf.create_pdfDoc(pdf_type='table')
pdf.pdf_insert_table(table_headers, addfinalLine=True)

parameter_dict = {'Redshift': z_objs[0],
                  'Programme ID': '102-B'}

for param, value in parameter_dict.items():

    row = ['-'] * len(table_headers)

    row[0] = param
    row[1] = value

    pdf.addTableRow(row, last_row=False)

pdf.table.add_hline()
pdf.generate_pdf(tablesFolder/'CGCG007_properties')