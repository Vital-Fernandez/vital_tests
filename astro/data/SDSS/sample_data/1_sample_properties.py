import os
import pyneb as pn
import numpy as np
import pandas as pd
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from pathlib import Path
from src.specsiser.data_printing import PdfPrinter


def list_files(directory, extension):
    output_list = []
    for file in os.listdir(directory):
        if file.endswith(extension):
            output_list.append(os.path.join(directory, file))
    return output_list

def colorChooser(ObsRatio, TheRatio):

    if ObsRatio > TheRatio:
        color = 'black'

    elif (TheRatio  < ObsRatio * 0.95):
        color = 'ForestGreen'  # 'green'#

    elif (TheRatio < ObsRatio * 0.85 ):
        color = 'YellowOrange'  # 'yellow'#

    else:
        color = 'YellowOrange'  # 'yellow'#

    return color


# Import the observation data
obsData = sr.loadConfData('D:/Pycharm Projects/vital_tests/astro/data/SDSS/flux_comparison.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['file_information']['data_folder'])
fileList = list_files(data_folder, '.fits')
addressList = list(data_folder / file for file in fileList)
flux_norm = obsData['sample_data']['norm_flux']
tableHeaders = [r'Object name', r'$EW(H\beta)$ $(\AA)$', r'$c(H\beta)$', r'$F(H\alpha)/F(H\beta)$']

scaleTable = 1000.0

fitsFolder, fitsFile = addressList[0].parent, addressList[0].name
sampleFolder, tableName = fitsFolder/'sample_data', 'sample_properties'

pdf = PdfPrinter()
pdf.create_pdfDoc(sampleFolder / tableName, pdf_type='table')
pdf.pdf_insert_table(tableHeaders)

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Open lineslog
    fitsFolder, fitsFile = file_address.parent, file_address.name
    lineLogFolder, lineLogFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_linesLog.txt')
    pdfLogFolder, pdfLogFile = fitsFolder / 'flux_analysis', fitsFile.replace('.fits', '_linesLog')

    objName = fitsFile.replace('.fits', '')
    print(f'- {i}: {objName}')

    # Set and crop the wavelength
    wave_rest, flux, header = sr.import_fits_data(fitsFolder/fitsFile, instrument='SDSS')
    idx_wave = (wave_rest >= obsData['sample_data']['wmin_array']) & (wave_rest <= obsData['sample_data']['wmax_array'])

    # Load line measurer object
    lm = sr.LineMesurerGUI(wave_rest[idx_wave], flux[idx_wave], lineLogFolder/lineLogFile, normFlux=flux_norm)

    # Measure line fluxes
    idcs_lines = ~lm.linesDF.index.str.contains('_b')
    obsLines = lm.linesDF.loc[idcs_lines].index.values

    # Equivalent width Hbeta
    eqw_Hbeta, eqwErr_Hbeta = lm.linesDF.loc['H1_4861A', 'eqw'], lm.linesDF.loc['H1_4861A', 'eqw_err']
    eqw_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(eqw_Hbeta, eqwErr_Hbeta)

    # Normalizing flux
    flux_Halpha = lm.linesDF.loc['H1_6563A', 'gauss_flux']
    flux_Hbeta = lm.linesDF.loc['H1_4861A', 'intg_flux']
    halpha_norm = flux_Halpha / flux_Hbeta
    halphaRatio_entry = r'${:0.2f}$'.format(halpha_norm)
    colorGrade = colorChooser(halpha_norm, 2.86)
    RatioComparison = r'\textcolor{{{}}}{{{}}}'.format(colorGrade, halphaRatio_entry)

    # Reddening coefficient
    rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
    rc.setCorr(obs_over_theo=halpha_norm/2.86, wave1=6563., wave2=4861.)
    cHbeta_entry = r'${:0.2f}$'.format(rc.cHbeta)
    pdf.addTableRow([objName, eqw_entry, cHbeta_entry, RatioComparison], last_row=True if file_address == addressList[-1] else False)

# Generate the table
pdf.generate_pdf(clean_tex=True)
