import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path

objList = ['gp030321', 'gp101157', 'gp121903']
conf_file_address = '../gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList=objList, group_variables=False)

dataFolder = Path(obsData['file_information']['data_folder'])
fileList = obsData['file_information']['files_list']

objList_B = obsData['file_information']['objectB_list']
fileList_B = obsData['file_information']['filesB_list']
objList_R = obsData['file_information']['objectR_list']
fileList_R = obsData['file_information']['filesR_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])

counter = 0
for i, obj in enumerate(objList):

    z = z_objs[i]
    fit_conf = obsData[f'{obj}_line_fitting']

    for ext in ('_BR', '_B', '_R'):

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objMask = dataFolder/'flux_analysis'/f'{obj}{ext}_mask.txt'
        lineLog_file = dataFolder/'flux_analysis'/f'{obj}{ext}_linesLog.txt'
        lineGrid_file = dataFolder/'flux_analysis'/f'{obj}{ext}_lineGrid.png'

        # Set wavelength and flux
        print(f'\n-- Treating {counter} :{obj}{ext}.fits')
        wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        wave_rest = wave / (1 + z)
        idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])
        flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array

        # Load line measurer object
        maskDF = pd.read_csv(objMask, delim_whitespace=True, header=0, index_col=0)
        lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave], normFlux=flux_norm)

        # Fit and check the regions
        obsLines = maskDF.index.values
        for j, lineLabel in enumerate(obsLines):

            print(f'-- {lineLabel}:')
            wave_regions = maskDF.loc[lineLabel, 'w1':'w6'].values
            lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=fit_conf)

            if lm.blended_check:
                plotFile = f'{obj}{ext}_{lineLabel}.png'
                lm.plot_fit_components(lm.fit_output, output_address=dataFolder/'flux_analysis'/plotFile)
        lm.save_lineslog(lm.linesDF, lineLog_file)

        # Plot the single lines:
        idcs_unblended = ~lm.linesDF.index.str.contains('_b')
        lm.plot_line_grid(lm.linesDF.loc[idcs_unblended], ncols=8, output_address=lineGrid_file)

        # Increase counter for obj number
        counter += 1
