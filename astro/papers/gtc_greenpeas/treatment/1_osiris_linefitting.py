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
reduction_dict = ['GRISM', 'FILTER4', 'OSIVERS', 'AIRMASS', 'EXPTIME']

ext = 'BR'
cycle = 'it1'

for i, obj in enumerate(objList):

    print(f'\n- Treating {i} :{obj}_{ext}.fits')

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    fits_file = dataFolder / f'{obj}_{ext}.fits'
    objMask = dataFolder / f'{obj}_{ext}_mask.txt'    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    fits_file = dataFolder / f'{obj}_{ext}.fits'
    objMask = dataFolder / f'{obj}_{ext}_mask.txt'

    # Declare output files
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
    lineGrid_file = objFolder / f'{obj}_{ext}_lineGrid_{cycle}.png'
    lineTable_file = objFolder / f'{obj}_{ext}_lineTable_{cycle}'

    # Load data
    wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
    flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
    maskDF = pd.read_csv(objMask, delim_whitespace=True, header=0, index_col=0)

    # Create line measurer object
    lm = sr.LineMesurer(wave, flux, redshift=z_objs[i], normFlux=flux_norm, crop_waves=(wmin_array[i], wmax_array[i]))
    lm.plot_spectrum_components()

    # # # Fit and check the regions
    # obsLines = maskDF.index.values
    # for j, lineLabel in enumerate(obsLines):
    #
    #     print(f'-- {lineLabel}:')
    #     wave_regions = maskDF.loc[lineLabel, 'w1':'w6'].values
    #     lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=obsData[f'{obj}_line_fitting'])
    #     # lm.print_results(show_fit_report=True, show_plot=True)
    #
    #     if lm.blended_check:
    #         plotFile = f'{obj}_{ext}_deblend_{lineLabel}_{cycle}.png'
    #         lm.plot_fit_components(lm.fit_output, output_address=objFolder/plotFile)
    #
    # # Save the lines log
    # lm.save_lineslog(lm.linesDF, lineLog_file)
    #
    # # Plot the single lines:
    # idcs_unblended = ~lm.linesDF.index.str.contains('_b')
    # lm.plot_line_grid(lm.linesDF.loc[idcs_unblended], ncols=8, output_address=lineGrid_file)
    #
    # # Lines to store in tables
    # idcs_lines = ~lm.linesDF.index.str.contains('_b')
    # linesLogOutput_df = lm.linesDF.loc[idcs_lines]
    #
    # # Reddening correction for flux tables
    # cHbeta, rc_pyneb = red_corr_HalphaHbeta_ratio(linesLogOutput_df, 0.0)
    #
    # # Table for the output data
    # print(f'- Printing results tables')
    # lm.table_fluxes(linesLogOutput_df, lineTable_file, rc_pyneb)
