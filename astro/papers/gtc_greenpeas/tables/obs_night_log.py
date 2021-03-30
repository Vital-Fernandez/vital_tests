import src.specsiser as sr
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from src.specsiser.data_printing import PdfPrinter
import pyneb as pn

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)
objList = obsData['file_information']['object_list']
objRef = obsData['file_information']['refName_list']
tables_folder = Path(obsData['file_information']['tables_folder'])

table_headers = {}

table_headers = ['SDSS ID', 'Label', 'R.A', 'DEC', 'z', '$Grisms^{1}$', '$T_{exp}$', 'Seeing', r'Standard stars$^{2}$']

pdf = PdfPrinter()
pdf.create_pdfDoc(pdf_type=None)
pdf.pdf_insert_table(table_headers, addfinalLine=True)

for i, obj in enumerate(objList):

    if i < 3:

        row = ['-'] * len(table_headers)
        objData = obsData[obj]

        # row[0] = objData['SDSS_label']
        row[0] = r'\href{{{}}}{{{}}}'.format(objData['SDSS_website'].replace('&', '\&'), objData['SDSS_label'])
        print(row[0])
        row[1] = objRef[i]
        row[2] = objData['RA']
        row[3] = objData['DEC']
        row[4] = objData['z']
        row[5] = f'{objData["grism_b"]}, {objData["grism_r"]}'
        # row[5] = f'{objData["grism_b"]} (blue), {objData["grism_r"]}(red)'
        row[6] = r'${}\times{}s$'.format(objData["exp_number"], objData["exp_time"])
        row[7] = r'$\approx {}$'.format(objData["seeing"])
        row[8] = f'{objData["standard_blue"]}, {objData["standard_red"]}'
        # row[8] = f'{objData["standard_blue"]} (blue), {objData["standard_red"]}(red)'
        # row[9] = f'{objData["night_type"]} night'

        pdf.addTableRow(row, last_row=False)

pdf.table.add_hline()
pdf.generate_pdf(tables_folder/'night_log')