import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import red_corr_HalphaHbeta_ratio

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
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

    for ext in ['_SDSS']:

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objFolder = resultsFolder/f'{obj}'
        objMask = objFolder/f'{obj}{ext}_mask.txt'
        lineLog_file, lineGrid_file = objFolder/f'{obj}{ext}_linesLog.txt', objFolder/f'{obj}{ext}_lineGrid.png'
        pdfTableFile, txtTableFile = objFolder/f'{obj}{ext}_linesTable', objFolder/f'{obj}{ext}_linesTable.txt'

        # if obj in ['gp004054', 'gp113303', 'gp232539']:
        if fits_file.is_file():

            # Load spectrum
            print(f'\n-- Treating {counter} :{obj}{ext}.fits')
            wave, flux_dict, header = sr.import_fits_data(fits_file, instrument='SDSS')
            flux, err = flux_dict['flux'], flux_dict['ivar']

            # Load line measurer object
            maskDF = pd.read_csv(objMask, delim_whitespace=True, header=0, index_col=0)
            lm = sr.LineMesurer(wave, flux, input_err=(1/err)/flux_norm, redshift=z, normFlux=flux_norm)
            # lm.plot_spectrum_components()

            # # Identify the emission lines
            # cont_norm_flux = lm.continuum_remover(noise_region)
            # obsLinesTable = lm.line_finder(cont_norm_flux, noiseWaveLim=noise_region, intLineThreshold=3)
            # obsLinesDF = lm.match_lines(obsLinesTable, maskDF, find_line_borders=False)
            # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=obsLinesDF)
            # # idcsObsLines = (obsLinesDF.observation == 'detected')
            # # lm.save_lineslog(obsLinesDF.loc[idcsObsLines], objFolder/objMask)

            # Fit and check the regions
            obsLines = maskDF.index.values
            fit_conf = obsData['SDSS_line_fitting']
            for j, lineLabel in enumerate(obsLines):

                print(f'-- {lineLabel}:')
                wave_regions = maskDF.loc[lineLabel, 'w1':'w6'].values
                lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=fit_conf)

                if lm.blended_check:
                    plotFile = f'{obj}{ext}_{lineLabel}.png'
                    # lm.print_results(show_fit_report=True, show_plot=True)
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
