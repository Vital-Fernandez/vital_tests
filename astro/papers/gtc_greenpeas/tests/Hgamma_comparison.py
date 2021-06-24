import numpy as np
import src.specsiser as sr
from pathlib import Path
import matplotlib.pyplot as plt
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
idx_band = int(obsData['file_information']['band_flux'])

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']

i_obj = 2
obj = objList[i_obj]
ext = '_BR'
cycle = 'c2'

for i_obj, obj in enumerate(objList):
    if i_obj < 3:
        print(f'Treating: {obj}')
        objFolder = resultsFolder / f'{obj}'
        lineLog_file = objFolder / f'{obj}{ext}_linesLog_{cycle}.txt'
        linesDF = sr.lineslogFile_to_DF(lineLog_file)

        fits_file = dataFolder / f'{obj}{ext}.fits'
        wave, flux, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        lm = sr.LineMesurer(wave, flux, redshift=z_objs[i_obj], normFlux=flux_norm, crop_waves=(wmin_array[i_obj], wmax_array[i_obj]))

        fits_blue_file = dataFolder / f'{obj}_B.fits'
        wave_b, flux, header_b = sr.import_fits_data(fits_blue_file, instrument='OSIRIS')
        flux_b = flux[idx_band][0]
        lm_b = sr.LineMesurer(wave_b, flux_b, redshift=z_objs[i_obj], normFlux=flux_norm, crop_waves=(wmin_array[i_obj], wmax_array[i_obj]))

        fits_red_file = dataFolder / f'{obj}_R.fits'
        wave_r, flux, header_r = sr.import_fits_data(fits_red_file, instrument='OSIRIS')
        flux_r = flux[idx_band][0]
        lm_r = sr.LineMesurer(wave_r, flux_r, redshift=z_objs[i_obj], normFlux=flux_norm, crop_waves=(wmin_array[i_obj], wmax_array[i_obj]))

        lineLabel = 'H1_6563A'
        lm.fit_from_wavelengths(lineLabel, line_wavelengths=linesDF.loc[lineLabel, 'w1':'w6'].values)
        # lm.print_results(show_plot=True)
        print(lineLabel, 'Joined', lm.intg_flux, lm.intg_err)

        lm_b.fit_from_wavelengths(lineLabel, line_wavelengths=linesDF.loc[lineLabel, 'w1':'w6'].values)
        # lm_b.print_results(show_plot=True)
        print(lineLabel, 'Blue', lm_b.intg_flux, lm_b.intg_err)

        print(f'Chang in {lineLabel} is b -> p: {1 - lm_b.intg_flux / lm.intg_flux}')
        # # Plot spectra components
        # fig, ax = plt.subplots(figsize=(12, 8))
        # ax.plot(lm.wave, lm.flux, label='Object combined spectrum', linestyle='--', color='tab:purple')
        # ax.plot(lm_b.wave, lm_b.flux, label='Object blue spectrum', color='tab:blue')
        # ax.plot(lm_r.wave, lm_r.flux, label='Object red spectrum', color='tab:red')
        # ax.legend()
        # ax.set_yscale('log')
        # plt.tight_layout()
        # plt.show()
