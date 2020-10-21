import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from fitelp.read_spectra import read_spectra
from matplotlib import pyplot as plt
from uncertainties.unumpy import uarray


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
        output_log = data_folder/f'{obj}_linesLog.txt'

        # Blue and red arm have a different data array structure
        if obj == 'B6479s':
            wave, flux = wave_data[0], flux_data[0][0]
        if obj == 'R8731s':
            wave, flux = wave_data[0], flux_data[0]

        # Individual line measurement
        lm = sr.LineMesurer(wave, flux, input_err=flux_data[0][1], normFlux=norm_flux, redshift=z_mean)
        fitConf = obsData[f'{obj}_individual_line_fitting']
        lm.plot_spectrum_components(continuumFlux=lm.errFlux)

        # Loop through the lines
        gauss_fluxes = {}
        sigma_fluxes = {}
        for lineLabel in ('H1_4861A_b', 'O3_4959A_b', 'O3_5007A_b'):

            lineWaves = mask_DF.loc[lineLabel[0:-2], 'w1':'w6'].values

            # Perform fit
            lm.fit_from_wavelengths(lineLabel, lineWaves, fit_conf=fitConf)

            # Display results
            print(lm)
            gauss_fluxes[lineLabel] = lm.lineGaussFlux
            gauss_fluxes[lineLabel+'Err'] = lm.lineGaussErr
            sigma_fluxes[lineLabel] = lm.p1[2]
            sigma_fluxes[lineLabel+'Err'] = lm.p1_Err[2]

            # Plot fit
            outputFile = data_folder/f'{obj}_{lineLabel}_multifitting.png'
            lm.plot_fit_components(lm.fit_output)#, output_address=outputFile)

        u_4861 = uarray(gauss_fluxes['H1_4861A_b'], gauss_fluxes['H1_4861A_bErr'])
        u_5007 = uarray(gauss_fluxes['O3_5007A_b'], gauss_fluxes['O3_5007A_bErr'])
        u_4959 = uarray(gauss_fluxes['O3_4959A_b'], gauss_fluxes['O3_4959A_bErr'])

        print(u_5007/u_4861)
        print(u_5007/u_4959)
        print(sigma_fluxes['O3_5007A_b'], sigma_fluxes['O3_5007A_bErr'])
        print(sigma_fluxes['O3_4959A_b'], sigma_fluxes['O3_4959A_bErr'])

        # Multi oxygen lines measurement same sigma
        # lm = sr.LineMesurer(wave, flux, normFlux=norm_flux, redshift=z_mean)
        # fitConf = obsData[f'{obj}_2multi_line_fitting']
        #
        # lineLabel = 'O3_5007A_b'
        # gauss_fluxes = {}
        # sigma_angstroms = {}
        #
        # O3_4959A_lineWaves = mask_DF.loc['O3_4959A', 'w1':'w6'].values
        # O3_5007A_lineWaves = mask_DF.loc['O3_5007A', 'w1':'w6'].values
        #
        # lineWaves = np.empty(6)
        # lineWaves[0:2], lineWaves[4:6] = O3_4959A_lineWaves[0:2], O3_5007A_lineWaves[4:6]
        # lineWaves[2], lineWaves[3] = O3_4959A_lineWaves[2], O3_5007A_lineWaves[3]
        #
        # # Perform fit
        # lm.fit_from_wavelengths(lineLabel, lineWaves, fit_conf=fitConf)
        #
        # # Display results
        # print(lm)
        # gauss_fluxes['O3_5007A_b'] = np.array([lm.fit_output.params['O3_5007A_n1_amplitude'].value,
        #                                     lm.fit_output.params['O3_5007A_n2_amplitude'].value,
        #                                     lm.fit_output.params['O3_5007A_w1_amplitude'].value])
        # gauss_fluxes['O3_4959A_b'] = np.array([lm.fit_output.params['O3_4959A_n1_amplitude'].value,
        #                                     lm.fit_output.params['O3_4959A_n2_amplitude'].value,
        #                                     lm.fit_output.params['O3_4959A_w1_amplitude'].value])
        # gauss_fluxes['O3_5007A_bErr'] = np.array([lm.fit_output.params['O3_5007A_n1_amplitude'].stderr,
        #                                     lm.fit_output.params['O3_5007A_n2_amplitude'].stderr,
        #                                     lm.fit_output.params['O3_5007A_w1_amplitude'].stderr])
        # gauss_fluxes['O3_4959A_bErr'] = np.array([lm.fit_output.params['O3_4959A_n1_amplitude'].stderr,
        #                                     lm.fit_output.params['O3_4959A_n2_amplitude'].stderr,
        #                                     lm.fit_output.params['O3_4959A_w1_amplitude'].stderr])
        #
        # sigma_angstroms['O3_5007A_b'] = np.array([lm.fit_output.params['O3_5007A_n1_sigma'].value,
        #                                     lm.fit_output.params['O3_5007A_n2_sigma'].value,
        #                                     lm.fit_output.params['O3_5007A_w1_sigma'].value])
        # sigma_angstroms['O3_5007A_bErr'] = np.array([lm.fit_output.params['O3_5007A_n1_sigma'].stderr,
        #                                     lm.fit_output.params['O3_5007A_n2_sigma'].stderr,
        #                                     lm.fit_output.params['O3_5007A_w1_sigma'].stderr])
        # sigma_angstroms['O3_4959A_b'] = np.array([lm.fit_output.params['O3_4959A_n1_sigma'].value,
        #                                     lm.fit_output.params['O3_4959A_n2_sigma'].value,
        #                                     lm.fit_output.params['O3_4959A_w1_sigma'].value])
        #
        # u_5007 = uarray(gauss_fluxes['O3_5007A_b'], gauss_fluxes['O3_5007A_bErr'])
        # u_4959 = uarray(gauss_fluxes['O3_4959A_b'], gauss_fluxes['O3_4959A_bErr'])
        # print(gauss_fluxes['O3_5007A_b'] / gauss_fluxes['O3_4959A_b'])
        # print(u_5007/u_4959)
        # print(sigma_angstroms['O3_5007A_b'], sigma_angstroms['O3_5007A_bErr'])
        # print(sigma_angstroms['O3_4959A_b'])
        #
        #
        # # Plot fit
        # outputFile = data_folder/f'{obj}_oxygens_multifitting.png'
        # lm.plot_fit_components(lm.fit_output)#, output_address=outputFile)


        # # Multi oxygen lines measurement theoretical flux
        # lm = sr.LineMesurer(wave, flux, normFlux=norm_flux, redshift=z_mean)
        # fitConf = obsData[f'{obj}_2multi_fixFlux_line_fitting']
        #
        # lineLabel = 'O3_5007A_b'
        # gauss_fluxes = {}
        # sigma_angstroms = {}
        #
        # O3_4959A_lineWaves = mask_DF.loc['O3_4959A', 'w1':'w6'].values
        # O3_5007A_lineWaves = mask_DF.loc['O3_5007A', 'w1':'w6'].values
        #
        # lineWaves = np.empty(6)
        # lineWaves[0:2], lineWaves[4:6] = O3_4959A_lineWaves[0:2], O3_5007A_lineWaves[4:6]
        # lineWaves[2], lineWaves[3] = O3_4959A_lineWaves[2], O3_5007A_lineWaves[3]
        #
        # # Perform fit
        # lm.fit_from_wavelengths(lineLabel, lineWaves, fit_conf=fitConf)
        #
        # # Display results
        # print(lm)
        # gauss_fluxes['O3_5007A_b'] = np.array([lm.fit_output.params['O3_5007A_n1_amplitude'].value,
        #                                     lm.fit_output.params['O3_5007A_n2_amplitude'].value,
        #                                     lm.fit_output.params['O3_5007A_w1_amplitude'].value])
        # gauss_fluxes['O3_4959A_b'] = np.array([lm.fit_output.params['O3_4959A_n1_amplitude'].value,
        #                                     lm.fit_output.params['O3_4959A_n2_amplitude'].value,
        #                                     lm.fit_output.params['O3_4959A_w1_amplitude'].value])
        # gauss_fluxes['O3_5007A_bErr'] = np.array([lm.fit_output.params['O3_5007A_n1_amplitude'].stderr,
        #                                     lm.fit_output.params['O3_5007A_n2_amplitude'].stderr,
        #                                     lm.fit_output.params['O3_5007A_w1_amplitude'].stderr])
        # gauss_fluxes['O3_4959A_bErr'] = np.array([lm.fit_output.params['O3_4959A_n1_amplitude'].stderr,
        #                                     lm.fit_output.params['O3_4959A_n2_amplitude'].stderr,
        #                                     lm.fit_output.params['O3_4959A_w1_amplitude'].stderr])
        #
        # print(gauss_fluxes['O3_5007A_b'] / gauss_fluxes['O3_4959A_b'])
        #
        # sigma_angstroms['O3_5007A_b'] = np.array([lm.fit_output.params['O3_5007A_n1_sigma'].value,
        #                                     lm.fit_output.params['O3_5007A_n2_sigma'].value,
        #                                     lm.fit_output.params['O3_5007A_w1_sigma'].value])
        # sigma_angstroms['O3_4959A_b'] = np.array([lm.fit_output.params['O3_4959A_n1_sigma'].value,
        #                                     lm.fit_output.params['O3_4959A_n2_sigma'].value,
        #                                     lm.fit_output.params['O3_4959A_w1_sigma'].value])
        #
        # sigma_angstroms['O3_5007A_bErr'] = np.array([lm.fit_output.params['O3_5007A_n1_sigma'].stderr,
        #                                     lm.fit_output.params['O3_5007A_n2_sigma'].stderr,
        #                                     lm.fit_output.params['O3_5007A_w1_sigma'].stderr])
        # sigma_angstroms['O3_4959A_bErr'] = np.array([lm.fit_output.params['O3_4959A_n1_sigma'].stderr,
        #                                     lm.fit_output.params['O3_4959A_n2_sigma'].stderr,
        #                                     lm.fit_output.params['O3_4959A_w1_sigma'].stderr])
        #
        # print(sigma_angstroms['O3_5007A_b'], sigma_angstroms['O3_5007A_bErr'])
        # print(sigma_angstroms['O3_4959A_b'], sigma_angstroms['O3_4959A_bErr'])
        #
        # # Plot fit
        # outputFile = data_folder/f'{obj}_oxygens_multifitting.png'
        # lm.plot_fit_components(lm.fit_output)#, output_address=outputFile)




        # Multi oxygen lines + Hbeta measurement
        # lm = sr.LineMesurer(wave, flux, normFlux=norm_flux, redshift=z_mean)
        # fitConf = obsData[f'{obj}_3multi_line_fitting']
        #
        # lineLabel = 'O3_5007A_b'
        #
        # H1_4861A_lineWaves = mask_DF.loc['H1_4861A', 'w1':'w6'].values
        # O3_5007A_lineWaves = mask_DF.loc['O3_5007A', 'w1':'w6'].values
        #
        # lineWaves = np.empty(6)
        # lineWaves[0:2], lineWaves[4:6] = H1_4861A_lineWaves[0:2], O3_5007A_lineWaves[4:6]
        # lineWaves[2], lineWaves[3] = H1_4861A_lineWaves[2], O3_5007A_lineWaves[3]
        #
        # # Perform fit
        # lm.fit_from_wavelengths(lineLabel, lineWaves, fit_conf=fitConf)
        #
        # # Display results
        # print(lm)
        #
        # # Plot fit
        # # outputFile = data_folder/f'{obj}_Hbeta-oxygens_multifitting.png'
        # lm.plot_fit_components(lm.fit_output)#, output_address=outputFile)
