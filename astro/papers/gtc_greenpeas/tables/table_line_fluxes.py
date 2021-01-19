import numpy as np
import pandas as pd
import pyneb as pn
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import red_corr_HalphaHbeta_ratio
from src.specsiser.data_printing import PdfPrinter, label_decomposition, format_for_table
from pylatex import Document, Figure, NewPage, NoEscape, Package, Tabular, Section, Tabu, Table, LongTable, MultiColumn, MultiRow
from functools import partial

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
tables_folder = Path(obsData['file_information']['tables_folder'])
fileList = obsData['file_information']['files_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])
counter = 0

# Start table
pdfTableFile, txtTableFile = tables_folder/f'emissionLinesTable', tables_folder/f'emissionLinesTable.txt'

# Loop through the objects and get the existing lines:

tableDF = pd.DataFrame(columns=['wavelength', 'latexLabel'])
df_dict = {}
Hbeta_dict = {}

for i, obj in enumerate(objList):

    for ext in ['_BR']:

        if i < 3:

            # Declare files location
            fits_file = dataFolder / f'{obj}{ext}.fits'
            objFolder = resultsFolder / f'{obj}'

            lineLog_file = objFolder / f'{obj}{ext}_linesLog.txt'
            linesDF_i = sr.lineslogFile_to_DF(lineLog_file)
            df_dict[obj] = linesDF_i
            Hbeta_dict[obj] = linesDF_i.loc['H1_4861A', 'intg_flux']

            # Add new lines to the master df
            for line in linesDF_i.index:
                if (line not in tableDF.index) and ('_b' not in line):
                    tableDF.loc[line] = [linesDF_i.loc[line, 'wavelength'], linesDF_i.loc[line, 'latexLabel']]

#
tableDF.sort_values('wavelength', inplace=True)

print(tableDF)

# Reddening law
rc = pn.RedCorr(R_V=obsData['sample_data']['RV'], law=obsData['sample_data']['red_law'])
X_Hbeta, Xx_ref, Xx = rc.X(4861.0), rc.X(4861.0), rc.X(tableDF.wavelength.values)
f_lines = Xx / Xx_ref - 1
f_Hbeta = X_Hbeta / Xx_ref - 1

pdf_geometry_options = {'right': '1cm', 'left': '1cm', 'top': '1cm', 'bottom': '2cm', 'landscape': 'true',
                        'paperwidth': '30in', 'paperheight': '30in'}

pdfDoc = Document(pdfTableFile, documentclass=u'article', geometry_options=pdf_geometry_options)
pdfDoc.packages.append(Package('preview', options=['active', 'tightpage', ]))
pdfDoc.packages.append(Package('hyperref', options=['unicode=true', ]))
pdfDoc.append(NoEscape(r'\pagenumbering{gobble}'))
pdfDoc.packages.append(Package('nicefrac'))
pdfDoc.packages.append(Package('color', options=['usenames', 'dvipsnames', ]))
pdfDoc.append(NoEscape(r'\begin{preview}'))

table1 = Tabular('lccccccc')

table1.add_hline()
row_headers = ['',
               '',
               MultiColumn(2, align='c', data=objList[0].upper()),
               MultiColumn(2, align='c', data=objList[1].upper()),
               MultiColumn(2, align='c', data=objList[2].upper())]
table1.add_row(row_headers, escape=False, strict=False)

row_subHeaders = ['Line label', r'$f_{lambda}$',
                  r'$F(\lambda)$', r'$I(\lambda)$',
                  r'$F(\lambda)$', r'$I(\lambda)$',
                  r'$F(\lambda)$', r'$I(\lambda)$']
table1.add_row(row_subHeaders, escape=False, strict=False)
table1.add_hline()

for i, label in enumerate(tableDF.index):

    row_data = [''] * 8

    row_data[0] = tableDF.loc[label, 'latexLabel']
    row_data[1] = f'{f_lines[i]:.2f}'
    for j, obj in enumerate(objList):
        if j < 3:
            if label in df_dict[obj].index:
                row_data[2 + 2*j] = df_dict[obj].loc[label].obsFlux
                row_data[2 + 2*j + 1] = df_dict[obj].loc[label].obsInt
            else:
                row_data[2 + 2*j] = '-'
                row_data[2 + 2*j + 1] = '-'

    output_row = list(map(partial(format_for_table, rounddig=3), row_data))
    table1.add_row(output_row, escape=False, strict=False)
table1.add_hline()

parameterers_row = [r'$c(H\beta)$',
                    '',
                    MultiColumn(2, align='c', data=0.0),
                    MultiColumn(2, align='c', data=0.0),
                    MultiColumn(2, align='c', data=0.0)]
table1.add_hline()
pdfDoc.append(table1)
pdfDoc.append(NoEscape(r'\end{preview}'))
print(pdfTableFile)
pdfDoc.generate_pdf(clean_tex=False)

# row_headers = [lineID,
# 	       f_lambda,
#                MultiColumn(2, align='c', data='GP01'),
#                MultiColumn(2, align='c', data='GP02'),
#                MultiColumn(2, align='c', data='GP03')]
# table1.add_row(row_headers)
# table1.add_hline()
#
# row_subheaders = [lineID, f_lambda,
#                   F_lambda, I_lambda,
#                   F_lambda, I_lambda,
# 		   F_lambda, I_lambda]
# table1.add_row(row_subheaders)
# table1.add_hline()
#
# row_data = [label, f_lambda_value, F_lambda_value, I_lambda_value, F_lambda_value, I_lambda_value, F_lambda_value, I_lambda_value]
# table1.add_row(row_data)
#
# table1.add_hline()
# parameterers_row = [Parameter,
#  		    '',
# 		    MultiColumn(2, align='c', data='Parameter_value_GP01'),
# 			MultiColumn(2, align='c', data='Parameter_value_GP02'),
# 			MultiColumn(2, align='c', data='Parameter_value_GP03')]
# table1.add_hline()
#
# pdfDoc.append(table1)
# pdfDoc.generate_pdf(clean_tex=False)

# table_format = 'l' + 'c' * (len(column_headers) - 1)
#
# # Initiate the table
# with self.pdfDoc.create(Tabu(table_format)) as self.table:
#     if column_headers != None:
#         self.table.add_hline()
#         self.table.add_row(list(map(str, column_headers)), escape=False, strict=False
#         self.table.add_hline()
# # flux_Hbeta = lines_df.loc['H1_4861A', 'intg_flux']
#
#
# obsLines = lines_df.index.values
# for lineLabel in obsLines:
#
#     label_entry = lines_df.loc[lineLabel, 'latexLabel']
#     wavelength = lines_df.loc[lineLabel, 'wavelength']
#     eqw, eqwErr = lines_df.loc[lineLabel, 'eqw'], lines_df.loc[lineLabel, 'eqw_err']
#
#     flux_intg = lines_df.loc[lineLabel, 'intg_flux'] / flux_Hbeta * scaleTable
#     flux_intgErr = lines_df.loc[lineLabel, 'intg_err'] / flux_Hbeta * scaleTable
#     flux_gauss = lines_df.loc[lineLabel, 'gauss_flux'] / flux_Hbeta * scaleTable
#     flux_gaussErr = lines_df.loc[lineLabel, 'gauss_err'] / flux_Hbeta * scaleTable
#
#     if (lines_df.loc[lineLabel, 'blended'] != 'None') and ('_m' not in lineLabel):
#         flux, fluxErr = flux_gauss, flux_gaussErr
#         label_entry = label_entry + '$_{gauss}$'
#     else:
#         flux, fluxErr = flux_intg, flux_intgErr
#
#     # Correct the flux
#     corr = pyneb_rc.getCorrHb(wavelength)
#     intensity, intensityErr = flux * corr, fluxErr * corr
#
#     eqw_entry = r'${:0.2f}\,\pm\,{:0.2f}$'.format(eqw, eqwErr)
#     flux_entry = r'${:0.2f}\,\pm\,{:0.2f}$'.format(flux, fluxErr)
#     intensity_entry = r'${:0.2f}\,\pm\,{:0.2f}$'.format(intensity, intensityErr)
#
#     # Add row of data
#     tex_row_i = [label_entry, eqw_entry, flux_entry, intensity_entry]
#     txt_row_i = [label_entry, eqw, eqwErr, flux, fluxErr, intensity, intensityErr]
#
#     lastRow_check = True if lineLabel == obsLines[-1] else False
#     pdf.addTableRow(tex_row_i, last_row=lastRow_check)
#     tableDF.loc[lineLabel] = txt_row_i[1:]
#
# # Data last rows
# row_Hbetaflux = [r'$H\beta$ $(erg\,cm^{-2} s^{-1} \AA^{-1})$',
#                  '',
#                  flux_Hbeta,
#                  flux_Hbeta * pyneb_rc.getCorr(4861)]
#
# row_cHbeta = [r'$c(H\beta)$',
#               '',
#               float(pyneb_rc.cHbeta),
#               '']
#
# pdf.addTableRow(row_Hbetaflux, last_row=False)
# pdf.addTableRow(row_cHbeta, last_row=False)
# tableDF.loc[row_Hbetaflux[0]] = row_Hbetaflux[1:] + [''] * 3
# tableDF.loc[row_cHbeta[0]] = row_cHbeta[1:] + [''] * 3
#
# # Format last rows
# pdf.table.add_hline()
# pdf.table.add_hline()
#
# # Save the pdf table
# try:
#     pdf.generate_pdf(clean_tex=True)
# except:
#     print('-- PDF compilation failure')
#
# # Save the txt table
# with open(txt_address, 'wb') as output_file:
#     string_DF = tableDF.to_string()
#     string_DF = string_DF.replace('$', '')
#     output_file.write(string_DF.encode('UTF-8'))
#
