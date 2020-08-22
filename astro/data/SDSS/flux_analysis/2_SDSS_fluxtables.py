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


# Import the observation data
obsData = sr.loadConfData('D:/Pycharm Projects/vital_tests/astro/data/SDSS/flux_comparison.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['file_information']['data_folder'])
fileList = list_files(data_folder, '.fits')
addressList = list(data_folder / file for file in fileList)
flux_norm = obsData['sample_data']['norm_flux']
tableHeaders = [r'$Transition$', '$EW(\AA)$', '$F(\lambda)$', '$I(\lambda)$']

scaleTable = 1000.0


# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Open lineslog
    fitsFolder, fitsFile = file_address.parent, file_address.name
    lineLogFolder, lineLogFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_linesLog.txt')
    pdfLogFolder, pdfLogFile = fitsFolder / 'flux_analysis', fitsFile.replace('.fits', '_linesLog')

    objName = fitsFile.replace('.fits', '')

    # Set and crop the wavelength
    wave_rest, flux, header = sr.import_fits_data(fitsFolder/fitsFile, instrument='SDSS')
    idx_wave = (wave_rest >= obsData['sample_data']['wmin_array']) & (wave_rest <= obsData['sample_data']['wmax_array'])

    # Load line measurer object
    lm = sr.LineMesurerGUI(wave_rest[idx_wave], flux[idx_wave], lineLogFolder/lineLogFile, normFlux=flux_norm)
    pdf = PdfPrinter()

    # Measure line fluxes
    idcs_lines = ~lm.linesDF.index.str.contains('_b')
    obsLines = lm.linesDF.loc[idcs_lines].index.values

    # Measure line fluxes
    pdf.create_pdfDoc(pdfLogFolder/pdfLogFile, pdf_type='table')
    pdf.pdf_insert_table(tableHeaders)

    # Normalizing flux
    flux_Halpha = lm.linesDF.loc['H1_6563A', 'intg_flux']
    flux_Hbeta = lm.linesDF.loc['H1_4861A', 'intg_flux']
    halpha_norm = flux_Halpha / flux_Hbeta

    rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
    rc.setCorr(obs_over_theo=halpha_norm/2.86, wave1=6563., wave2=4861.)
    print(f'-- cHbeta {rc.cHbeta}')

    for lineLabel in obsLines:

        label_entry = sr._linesDb.loc[lineLabel, 'latexLabel']
        wavelength = lm.linesDF.loc[lineLabel, 'wavelength']

        eqw, eqwErr = lm.linesDF.loc[lineLabel, 'eqw'], lm.linesDF.loc[lineLabel, 'eqw_err']

        flux_intg = lm.linesDF.loc[lineLabel, 'intg_flux'] / flux_Hbeta * scaleTable
        flux_intgErr = lm.linesDF.loc[lineLabel, 'intg_err'] / flux_Hbeta * scaleTable
        flux_gauss = lm.linesDF.loc[lineLabel, 'gauss_flux'] / flux_Hbeta * scaleTable
        flux_gaussErr = lm.linesDF.loc[lineLabel, 'gauss_err'] / flux_Hbeta * scaleTable

        if lm.linesDF.loc[lineLabel, 'blended'] != 'None':
            flux, fluxErr = flux_gauss, flux_gaussErr
            label_entry = label_entry + '$_{gauss}$'
        else:
            flux, fluxErr = flux_intg, flux_intgErr

        # Correct the flux
        corr = rc.getCorrHb(wavelength)
        intensity, intensityErr = flux * corr, fluxErr * corr

        # Format the entries
        eqw_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(eqw, eqwErr)
        flux_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(flux, fluxErr)
        intensity_entry = r'${:0.2f}$ $\pm$ ${:0.2f}$'.format(intensity, intensityErr)

        lastRow_check = True if lineLabel == obsLines[-1] else False
        pdf.addTableRow([label_entry, eqw_entry, flux_entry, intensity_entry], last_row=lastRow_check)

    # Data last rows
    row_Hbetaflux = [r'$H\beta$ $(erg\,cm^{-2} s^{-1} \AA^{-1})$', '', flux_Hbeta*flux_norm, flux_Hbeta*flux_norm*rc.getCorr(4861)]
    row_cHbeta = [r'$c(H\beta)$', '', rc.cHbeta, '']

    # Format last rows
    pdf.addTableRow(row_Hbetaflux, last_row=False)
    pdf.addTableRow(row_cHbeta, last_row=False)
    pdf.table.add_hline()
    pdf.table.add_hline()

    # Generate the table
    pdf.generate_pdf(clean_tex=True)

