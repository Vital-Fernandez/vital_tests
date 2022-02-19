from pathlib import Path
import src.specsiser as sr
from delete.data_printing import PdfPrinter
from pylatex import NoEscape
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)
objList = obsData['file_information']['object_list']
objRef = obsData['file_information']['refName_list']
tables_folder = Path('/home/vital/Dropbox/Astrophysics/Seminars/UniVapo 2021/')

table_headers = {}

table_headers = ['SDSS ID', 'Label', 'R.A', 'DEC', 'z', '$Grisms^{1}$', '$T_{exp}$', 'Seeing', r'Standard stars$^{2}$']

pdf = PdfPrinter()

pdf.create_pdfDoc(pdf_type='table')
pdf.pdfDoc.append(NoEscape('\definecolor{background}{rgb}{0.169, 0.169, 0.169}'))
pdf.pdfDoc.append(NoEscape('\definecolor{foreground}{rgb}{0.702, 0.780, 0.847}'))
pdf.pdfDoc.append(NoEscape(r'\arrayrulecolor{foreground}'))

pdf.pdf_insert_table(table_headers, addfinalLine=True, color_font='foreground')

for i, obj in enumerate(objList):

    if i < 3:
        row = ['-'] * len(table_headers)
        objData = obsData[obj]

        # row[0] = objData['SDSS_label']
        row[0] = r'\href{{{}}}{{{}}}'.format(objData['SDSS_website'].replace('&', '\&'), objData['SDSS_label'])
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

        # pdf.addTableRow(NoEscape(r'\rowcolor{background}'), last_row=False)
        pdf.addTableRow(row, last_row=False, color_font='foreground')

pdf.table.add_hline()
pdf.generate_pdf(tables_folder/'night_log_dark', clean_tex=False)
