import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import red_corr_HalphaHbeta_ratio
import matplotlib.pyplot as plt

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


ext = '_BR'
cycle = 'c2'
cycle_ref = 'First_cycle'

counter = 0
for i, obj in enumerate(objList):

    z = z_objs[i]
    wmin, wmax = wmin_array[i], wmax_array[i]
    fit_conf = obsData[f'{obj}_line_fitting']

    for ext in ['_BR']:

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objFolder = resultsFolder/f'{obj}'
        objMask = objFolder/f'{obj}{ext}_mask.txt'
        nebCompFile = objFolder / f'{obj}{ext}_NebFlux_{cycle}.txt'
        stellarFluxFile = objFolder / f'{obj}{ext}_stellarFlux_{cycle}.txt'
        lineLog_file = objFolder/f'{obj}{ext}_linesLog_{cycle}.txt'
        lineGrid_file = objFolder/f'{obj}{ext}_lineGrid_{cycle}.png'

        # Load spectrum
        print(f'\n- Treating {counter} :{obj}{ext}.fits')
        wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array

        # Remove the stellar component
        wave_neb, flux_nebular = np.loadtxt(nebCompFile, unpack=True)
        wave_star, flux_stellar = np.loadtxt(stellarFluxFile, unpack=True)
        flux_neb, flux_star = flux_nebular/flux_norm, flux_stellar/flux_norm

        # Load line measurer object
        lm = sr.LineMesurer(wave, flux, redshift=z, normFlux=flux_norm, crop_waves=(wmin, wmax))

        # fig, ax = plt.subplots(figsize=(12, 8))
        # # ax.plot(lm.wave, lm.flux, label='Object flux')
        # # # ax.plot(stellar_Wave, obj_input_flux, label='Input starlight spectrum', linestyle='--')
        # # ax.plot(wave_neb, flux_neb, label='Nebular flux')
        # # ax.plot(wave_star, flux_star, label='Stellar flux')
        # # ax.plot(wave_star, flux_star + flux_neb, label='Combined continuum', linestyle=':')
        # ax.plot(lm.wave, lm.flux, label='Object flux')
        # ax.plot(lm.wave, lm.flux - flux_star, label='No stellar', linestyle='--')
        #
        #
        # ax.legend()
        # # ax.set_yscale('log')
        # plt.tight_layout()
        # plt.show()

        # Remove the stellar component
        lm.flux = lm.flux - flux_star

        # Load lines mask
        maskDF = pd.read_csv(objMask, delim_whitespace=True, header=0, index_col=0)
        # lm.plot_spectrum_components(flux_neb/flux_norm, log_scale=True)

        # # Fit and check the regions
        obsLines = maskDF.index.values
        for j, lineLabel in enumerate(obsLines):

            print(f'-- {lineLabel}:')
            wave_regions = maskDF.loc[lineLabel, 'w1':'w6'].values
            lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=fit_conf)
            # lm.print_results(show_fit_report=True, show_plot=True)

            if lm.blended_check:
                plotFile = f'{obj}{ext}_{lineLabel}_{cycle}.png'
                lm.plot_fit_components(lm.fit_output, output_address=objFolder/plotFile)
                # if i > 1:
                #     lm.print_results(show_fit_report=True, show_plot=True)

        # Save the lines log
        lm.save_lineslog(lm.linesDF, lineLog_file)

        # Plot the single lines:
        idcs_unblended = ~lm.linesDF.index.str.contains('_b')
        lm.plot_line_grid(lm.linesDF.loc[idcs_unblended], ncols=8, output_address=lineGrid_file, log_check=False)
        # lm.plot_line_grid(lm.linesDF.loc[idcs_unblended], ncols=8, log_check=False)

        # Lines to store in tables
        # idcs_lines = ~lm.linesDF.index.str.contains('_b')
        # linesLogOutput_df = lm.linesDF.loc[idcs_lines]
        #
        # # Reddening correction for flux tables
        # cHbeta, rc_pyneb = red_corr_HalphaHbeta_ratio(linesLogOutput_df, 0.0)
        #
        # # Table for the output data
        # print(f'- Printing results tables')
        # lm.table_fluxes(linesLogOutput_df, pdfTableFile, txtTableFile, rc_pyneb)

