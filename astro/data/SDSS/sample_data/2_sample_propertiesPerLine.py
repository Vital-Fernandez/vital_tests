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

scaleTable = 1000.0

fitsFolder, fitsFile = addressList[0].parent, addressList[0].name

df_dict = {}
base_DF = None

# Load the first object
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
    lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave], lineLogFolder / lineLogFile, normFlux=flux_norm)

    # Measure line fluxes
    if i == 0:
        base_DF = lm.linesDF.copy()
        df_dict[objName] = lm.linesDF.copy()
    else:
        df_dict[objName] = lm.linesDF.copy()
        idcs_missingLines = ~lm.linesDF.index.isin(base_DF.index)
        missing_lines = lm.linesDF.loc[idcs_missingLines].index.values

        for line in missing_lines:
            base_DF.loc[line] = lm.linesDF.loc[line]


#Sort the data
objectList = list(df_dict.keys())
objectList.sort()

base_DF.sort_values(by=['wavelength'], inplace=True)
idcsLines = ~base_DF.index.str.contains('_b')
lineList = base_DF.loc[idcsLines].index.values

# Print the data
for lineLabel in lineList:

    sampleFolder, tableName = fitsFolder/'sample_data'/'properties_by_line', f'{lineLabel}_fluxes'

    pdf = PdfPrinter()
    pdf.create_pdfDoc(sampleFolder / tableName, pdf_type='table')
    refarenceTable = base_DF.loc[lineLabel, 'latexLabel']
    print(refarenceTable)

    tableHeaders = [refarenceTable, '$EW(\AA)$', '$F(\lambda)$', '$I(\lambda)$', r'$c(H\beta)$']
    pdf.pdf_insert_table(tableHeaders)

    for obj in objectList:

        objDF = df_dict[obj]

        if lineLabel in objDF.index:

            # Normalizing flux
            flux_Halpha = objDF.loc['H1_6563A', 'gauss_flux']
            flux_Hbeta = objDF.loc['H1_4861A', 'intg_flux']
            halpha_norm = flux_Halpha / flux_Hbeta

            rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
            rc.setCorr(obs_over_theo=halpha_norm / 2.86, wave1=6563., wave2=4861.)
            cHbeta_entry = '${:0.2f}$'.format(rc.cHbeta)

            # Eqw
            eqw, eqwErr = objDF.loc['H1_4861A', 'eqw'], objDF.loc['H1_4861A', 'eqw_err']
            eqw_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(eqw, eqwErr)

            flux_intg = objDF.loc[lineLabel, 'intg_flux'] / flux_Hbeta * scaleTable
            flux_intgErr = objDF.loc[lineLabel, 'intg_err'] / flux_Hbeta * scaleTable
            flux_gauss = objDF.loc[lineLabel, 'gauss_flux'] / flux_Hbeta * scaleTable
            flux_gaussErr = objDF.loc[lineLabel, 'gauss_err'] / flux_Hbeta * scaleTable

            if objDF.loc[lineLabel, 'blended_label'] != 'None':
                flux, fluxErr = flux_gauss, flux_gaussErr
            else:
                flux, fluxErr = flux_intg, flux_intgErr

            # Correct the flux
            wavelength = objDF.loc[lineLabel, 'wavelength']
            corr = rc.getCorrHb(wavelength)
            intensity, intensityErr = flux * corr, fluxErr * corr

            flux_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(flux, fluxErr)
            intensity_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(intensity, intensityErr)

        else:
            eqw_entry, flux_entry, intensity_entry, cHbeta_entry = '-', '-', '-', '-'

        lastRow_check = True if obj == objectList[-1] else False
        pdf.addTableRow([obj, eqw_entry, flux_entry, intensity_entry, cHbeta_entry], last_row=lastRow_check)

    # Generate the table
    pdf.generate_pdf(clean_tex=True)

