import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import red_corr_HalphaHbeta_ratio

objList = ['gp030321', 'gp101157', 'gp121903']
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList=objList, group_variables=False)

dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])

counter = 0
for i, obj in enumerate(objList):

    z = z_objs[i]
    wmin, wmax = wmin_array[i], wmax_array[i]
    fit_conf = obsData[f'{obj}_line_fitting']

    # for ext in ('_BR', '_B', '_R'):
    for ext in ('_BR', '_B', '_R'):

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objFolder = resultsFolder/f'{obj}'
        objMask = objFolder/f'{obj}{ext}_mask.txt'
        lineLog_file, lineGrid_file = objFolder/f'{obj}{ext}_linesLog.txt', objFolder/f'{obj}{ext}_lineGrid.png'
        pdfTableFile, txtTableFile = objFolder/f'{obj}{ext}_linesTable', objFolder/f'{obj}{ext}_linesTable.txt'

        # Load spectrum
        print(f'\n-- Treating {counter} :{obj}{ext}.fits')
        wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
        print(header['GRISM'])

        # Load line measurer object
        maskDF = pd.read_csv(objMask, delim_whitespace=True, header=0, index_col=0)
        lm = sr.LineMesurer(wave, flux, input_err=flux_array[3][0], redshift=z, normFlux=flux_norm, crop_waves=(wmin, wmax))
        lm.plot_spectrum_components(continuumFlux=lm.errFlux)

        # Fit and check the regions
        obsLines = maskDF.index.values
        for j, lineLabel in enumerate(obsLines):

            print(f'-- {lineLabel}:')
            wave_regions = maskDF.loc[lineLabel, 'w1':'w6'].values
            lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=fit_conf)
            print(lm)
            lm.plot_fit_components(lm.fit_output)

            if lm.blended_check:
                plotFile = f'{obj}{ext}_{lineLabel}.png'
                lm.plot_fit_components(lm.fit_output, output_address=objFolder/plotFile)

        # Save the lines log
        lm.save_lineslog(lm.linesDF, lineLog_file)

        # Plot the single lines:
        idcs_unblended = ~lm.linesDF.index.str.contains('_b')
        lm.plot_line_grid(lm.linesDF.loc[idcs_unblended], ncols=8, output_address=lineGrid_file)

        # Lines to store in tables
        idcs_lines = ~lm.linesDF.index.str.contains('_b')
        linesLogOutput_df = lm.linesDF.loc[idcs_lines]

        # Reddening correction for flux tables
        cHbeta, rc_pyneb = red_corr_HalphaHbeta_ratio(linesLogOutput_df, obsData[obj]['cHbeta'])

        # Table for the output data
        print(f'- Printing results tables')
        lm.table_fluxes(linesLogOutput_df, pdfTableFile, txtTableFile, rc_pyneb)

        # Increase counter for obj number
        counter += 1
