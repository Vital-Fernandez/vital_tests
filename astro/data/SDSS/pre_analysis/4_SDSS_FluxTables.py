import os
import numpy as np
import pandas as pd
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from pathlib import Path
from delete.data_printing import PdfPrinter

def list_files(directory, extension):
    output_list = []
    for file in os.listdir(directory):
        if file.endswith(extension):
            output_list.append(os.path.join(directory, file))
    return output_list


def colorChooser(ObsRatio, TheRatio):
    if (TheRatio * 0.95 < ObsRatio < TheRatio * 1.05):
        color = 'ForestGreen'  # 'green'#

    elif (TheRatio * 0.90 < ObsRatio < TheRatio * 1.10):
        color = 'YellowOrange'  # 'yellow'#

    else:
        color = 'BrickRed'

    return color


# Import the observation data
obsData = sr.loadConfData('../flux_comparison.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/database/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['file_information']['data_folder'])
fileList = list_files(data_folder, '.fits')
addressList = list(data_folder/file for file in fileList)
fluxNorm = obsData['sample_data']['norm_flux']
tableHeaders = [r'$\lambda(\AA)$', '$EW(\AA)$', '$F_{intg}(\AA)$', '$F_{gauss}(\AA)$']

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    print(f'- Treating: {file_address}')

    # Open lineslog
    linesLogAddress = str(file_address).replace('.fits', '_treatedlinesLog.txt')
    tableLogAddress = str(file_address).replace('.fits', '_treatedlinesLog')
    figureLogAddress = str(file_address).replace('.fits', '_treatedlinesLog.png')

    # Set and crop the wavelength
    wave_rest, flux, header = sr.import_fits_data(file_address, instrument='SDSS')
    idx_wave = (wave_rest >= obsData['sample_data']['wmin_array']) & (wave_rest <= obsData['sample_data']['wmax_array'])

    # Load line measurer object
    lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave] / fluxNorm, linesDF_address=linesLogAddress)
    pdf = PdfPrinter()

    # Measure line fluxes
    idcs_lines = ~lm.linesDF.index.str.contains('_b')
    obsLines = lm.linesDF.loc[idcs_lines].index.values

    pdf.create_pdfDoc(tableLogAddress, pdf_type='table')
    pdf.pdf_insert_table(tableHeaders)

    flux_Hbeta = lm.linesDF.loc['H1_4861A', 'intg_flux']

    for lineLabel in obsLines:

        label_entry = lm.linesDF.loc[lineLabel, 'latexLabel']
        eqw, eqwErr = lm.linesDF.loc[lineLabel, 'eqw'], lm.linesDF.loc[lineLabel, 'eqw_err']
        flux_intg, flux_intgErr = lm.linesDF.loc[lineLabel, 'intg_flux'], lm.linesDF.loc[lineLabel, 'intg_err']
        flux_gauss, flux_gaussErr = lm.linesDF.loc[lineLabel, 'gauss_flux'], lm.linesDF.loc[lineLabel, 'gauss_err']

        eqw_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(eqw, eqwErr)
        intgF_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(flux_intg/flux_Hbeta, flux_intgErr/flux_Hbeta)
        gaussF_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(flux_gauss/flux_Hbeta, flux_gaussErr/flux_Hbeta)

        if lm.linesDF.loc[lineLabel, 'blended_label'] == 'None':
            colorGrade = colorChooser(flux_intg/flux_gauss, 1.0)
        else:
            colorGrade = 'black'

        gaussF_color_entry = r'\textcolor{{{}}}{{{}}}'.format(colorGrade, gaussF_entry)
        lastRow_check = True if lineLabel == obsLines[-1] else False
        pdf.addTableRow([label_entry, eqw_entry, intgF_entry, gaussF_color_entry], last_row=lastRow_check)

    # Plot the matched lines:
    lm.plot_line_mask_selection(lm.linesDF, ncols=5, output_address=figureLogAddress)

    # Generate the table
    pdf.generate_pdf(clean_tex=True)

