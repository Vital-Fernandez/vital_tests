import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from fitelp.read_spectra import read_spectra
from matplotlib import pyplot as plt

obj_list = ['B6479s', 'R8731s']
data_folder = Path('D:/Google drive/Astrophysics/Datos/broad_multiComponent')
obsData = sr.loadConfData('broad_conf.ini', group_variables=False)
z_mean, z_err = sr.redshift_calculation(obsData['sample_data']['obs_waves'], obsData['sample_data']['emis_waves'])
norm_flux = obsData['sample_data']['norm_flux']

for idx_obj, obj in enumerate(obj_list):

    if obj == 'B6479s':

        fits_address = data_folder/f'{obj}.fits'
        mask_address = data_folder/f'{obj}_mask.txt'
        wave_data, flux_data = read_spectra(fits_address, scaleFlux=1)
        mask_DF = sr.lineslogFile_to_DF(mask_address)

        # Blue and red arm have a different data array structure
        if obj == 'B6479s':
            wave, flux = wave_data[0], flux_data[0][0]
        if obj == 'R8731s':
            wave, flux = wave_data[0], flux_data[0]

        # Line measuring object
        lm = sr.LineMesurer(wave, flux, normFlux=norm_flux, redshift=z_mean)

        # Declare line and line region
        lineLabel = 'O3_5007A_b'
        lineWaves = np.array([4966.75, 4982.46, 4997.58, 5013.89, 5022.57, 5040.09])

        # Fit configuration
        fit_conf = {'O3_5007A_b': 'O3_5007A_n1-O3_5007A_n2-O3_5007A_w1',
                    'O3_5007A_w1_sigma': {'value': 5.0, 'min': 0.0}}

        # Perform fit
        lm.fit_from_wavelengths(lineLabel, lineWaves, fit_conf=fit_conf)

        # Display results
        print(lm)

        # Plot fit
        lm.plot_fit_components(lm.fit_output)

        # # Line measuring object
        # lm = sr.LineMesurer(wave, flux, normFlux=norm_flux, redshift=z_mean)
        #
        # # Convert spectrum to velocity units
        # x_vel = ((lm.wave-5007)/5007) * 300000
        # lm.wave = x_vel
        #
        # # Declare line and line region
        # lineLabel = 'O3_0A_b'
        # lineWaves_vel = ((lineWaves-5007)/5007) * 300000
        #
        # # Fit configuration
        # fit_conf = {'O3_0A_b': 'O3_0A_n1-O3_0A_n2-O3_0A_w1',
        #
        #             'O3_0A_n1_center': {'value': 0.0},
        #             'O3_0A_n2_center': {'value': 0.0},
        #             'O3_0A_w1_center': {'value': 0.0},
        #
        #             'O3_0A_n1_sigma': {'value': 50.0, 'min': 0.0},
        #             'O3_0A_n2_sigma': {'value': 50.0, 'min': 0.0},
        #             'O3_0A_w1_sigma': {'value': 200.0, 'min': 100.0}}
        #
        # # Perform fit
        # lm.fit_from_wavelengths(lineLabel, lineWaves, fit_conf=fit_conf)
        #
        # # Display results
        # print(lm)
        #
        # # Plot fit
        # lm.plot_fit_components(lm.fit_output)

