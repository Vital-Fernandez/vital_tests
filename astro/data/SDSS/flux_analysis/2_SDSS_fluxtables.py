import os
import pyneb as pn
import numpy as np
import pandas as pd
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from pathlib import Path
from src.specsiser.data_printing import PdfPrinter
from astro.data.SDSS.shared_scripts import list_objName, obsConfaddress, obsFolder


# Declare data and files location
objList = list_objName(obsFolder, '.fits')
obsData = sr.loadConfData(obsConfaddress, objList, group_variables=False)
data_folder = Path(obsData['file_information']['data_folder'])
addressList = list(Path(f'{data_folder/objName}.fits') for objName in objList)

# Sample properties
norm_Flux = obsData['sample_data']['norm_flux']
wmin, wmax = obsData['sample_data']['wmin_array'], obsData['sample_data']['wmax_array']

# Table settings
tableHeaders = [r'$Transition$', '$EW(\AA)$', '$F(\lambda)$', '$I(\lambda)$']
scaleTable = 1000.0
dict_failures = {}

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Open lineslog
    objName = objList[i]
    fitsFolder, fitsFile = file_address.parent, file_address.name
    lineLogFolder, lineLogFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_linesLog.txt')
    pdfTableFolder, pdfTableFile = fitsFolder / 'flux_analysis', fitsFile.replace('.fits', '_linesTable')
    txtTableFolder, txtTableFile = fitsFolder / 'flux_analysis', fitsFile.replace('.fits', '_linesTable.txt')
    print(f'\n- {i}: {objName}')

    # Set and crop the wavelength
    wave_rest, flux, header = sr.import_fits_data(fitsFolder/fitsFile, instrument='SDSS')
    idx_wave = (wave_rest >= wmin) & (wave_rest <= wmax)

    # Load line measurer object
    lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave], lineLogFolder/lineLogFile, normFlux=norm_Flux)
    pdf = PdfPrinter()
    tableDF = pd.DataFrame(columns=tableHeaders[1:])

    # Measure line fluxes
    idcs_lines = ~lm.linesDF.index.str.contains('_b')
    obsLines = lm.linesDF.loc[idcs_lines].index.values

    # Measure line fluxes
    pdf.create_pdfDoc(pdfTableFolder/pdfTableFile, pdf_type='table')
    pdf.pdf_insert_table(tableHeaders)

    # Normalizing flux
    flux_Halpha = lm.linesDF.loc['H1_6563A', 'gauss_flux']
    flux_Hbeta = lm.linesDF.loc['H1_4861A', 'intg_flux']
    halpha_norm = flux_Halpha / flux_Hbeta

    rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
    rc.setCorr(obs_over_theo=halpha_norm/2.86, wave1=6563., wave2=4861.)

    for lineLabel in obsLines:

        label_entry = lm.linesDF.loc[lineLabel, 'latexLabel']
        wavelength = lm.linesDF.loc[lineLabel, 'wavelength']
        eqw, eqwErr = lm.linesDF.loc[lineLabel, 'eqw'], lm.linesDF.loc[lineLabel, 'eqw_err']

        flux_intg = lm.linesDF.loc[lineLabel, 'intg_flux'] / flux_Hbeta * scaleTable
        flux_intgErr = lm.linesDF.loc[lineLabel, 'intg_err'] / flux_Hbeta * scaleTable
        flux_gauss = lm.linesDF.loc[lineLabel, 'gauss_flux'] / flux_Hbeta * scaleTable
        flux_gaussErr = lm.linesDF.loc[lineLabel, 'gauss_err'] / flux_Hbeta * scaleTable

        if (lm.linesDF.loc[lineLabel, 'blended_label'] != 'None') and ('_m' not in lineLabel):
            flux, fluxErr = flux_gauss, flux_gaussErr
            label_entry = label_entry + '$_{gauss}$'
        else:
            flux, fluxErr = flux_intg, flux_intgErr

        # Correct the flux
        corr = rc.getCorrHb(wavelength)
        intensity, intensityErr = flux * corr, fluxErr * corr

        # Format the entries
        eqw_entry = r'${:0.2f}\pm{:0.2f}$'.format(eqw, eqwErr)
        flux_entry = r'${:0.2f}\pm{:0.2f}$'.format(flux, fluxErr)
        intensity_entry = r'${:0.2f}\pm{:0.2f}$'.format(intensity, intensityErr)

        # Add row of data
        row_i = [label_entry, eqw_entry, flux_entry, intensity_entry]
        lastRow_check = True if lineLabel == obsLines[-1] else False

        pdf.addTableRow([label_entry, eqw_entry, flux_entry, intensity_entry], last_row=lastRow_check)
        tableDF.loc[lineLabel] = row_i[1:]

    # Data last rows
    row_Hbetaflux = [r'$H\beta$ $(erg\,cm^{-2} s^{-1} \AA^{-1})$', '', flux_Hbeta*norm_Flux, flux_Hbeta*norm_Flux*rc.getCorr(4861)]
    row_cHbeta = [r'$c(H\beta)$', '', rc.cHbeta, '']

    pdf.addTableRow(row_Hbetaflux, last_row=False)
    pdf.addTableRow(row_cHbeta, last_row=False)
    tableDF.loc[row_Hbetaflux[0]] = row_Hbetaflux[1:]
    tableDF.loc[row_cHbeta[0]] = row_cHbeta[1:]

    # Format last rows
    pdf.table.add_hline()
    pdf.table.add_hline()

    # Save the pdf table
    try:
        pdf.generate_pdf(clean_tex=True)
    except:
        print('-- PDF compilation failure')
        dict_failures[i] = objName

    # Save the txt table
    with open(txtTableFolder/txtTableFile, 'wb') as output_file:
        string_DF = tableDF.to_string()
        string_DF = string_DF.replace('$', '')
        output_file.write(string_DF.encode('UTF-8'))

print('Total errors', len(dict_failures.keys()))
print(dict_failures.keys())
print(dict_failures.values())

