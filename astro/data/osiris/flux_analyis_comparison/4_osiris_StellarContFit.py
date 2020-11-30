import numpy as np
import pandas as pd
from pathlib import Path
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.starContinuum_functions import SSPsynthesizer, computeSSP_galaxy_mass
import matplotlib.pyplot as plt

objList = ['gp030321', 'gp101157', 'gp121903']
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList=objList, group_variables=False)
starlightFolder = Path(obsData['SSP_synthesis']['starlight_folder'])

fileList = obsData['file_information']['files_list']
dataFolder = Path(obsData['file_information']['data_folder'])
outputFolder = dataFolder/'flux_analysis'

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
        cHbeta = obsData[obj]['cHbeta']
        Te_low = obsData[obj]['Te_low']
        HeII_HII, HeIII_HeII = obsData[obj]['HeII_HII'], obsData[obj]['HeII_HII']

        for ext in ('_BR', '_B'):
        # for ext in ('_B'):

            # Declare files location
            fits_file = dataFolder/f'{obj}{ext}.fits'
            objMask = dataFolder/'flux_analysis'/f'{obj}{ext}_mask.txt'
            lineLog_file = outputFolder/f'{obj}{ext}_linesLog.txt'
            nebFluxNoNebCompFile = outputFolder/f'{obj}{ext}_obs_RemoveNebularComp.txt'
            massFracPlotFile = outputFolder/f'{obj}{ext}_SSP_MasFrac.png'
            LightFracPlotFile = outputFolder/f'{obj}{ext}_SSP_LightFrac.png'
            stellarPlotFile = outputFolder/f'{obj}{ext}_stellarFit.png'

            # Set wavelength and flux
            print(f'\n-- Treating {counter} :{obj}{ext}.fits')
            wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
            wave_rest = wave / (1 + z)
            idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])
            flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array

            # Load the data
            obsWave, obsFlux = wave_rest[idx_wave], flux[idx_wave]
            specWave, specFlux = np.loadtxt(nebFluxNoNebCompFile, unpack=True)

            # Measuring objects
            lm = sr.LineMesurer(specWave, specFlux, lineLog_file, normFlux=flux_norm)
            sw = SSPsynthesizer()

            # Generate starlight files
            idcs_lines = ~lm.linesDF.index.str.contains('_b')
            gridFileName, outputFile, saveFolder, waveResample, fluxResample = sw.generate_starlight_files(starlightFolder,
                                                                                                             f'{obj}{ext}',
                                                                                                             specWave,
                                                                                                             specFlux,
                                                                                                             lm.linesDF.loc[idcs_lines])

            # # Launch starlight
            # print(f'\n-Initiating starlight: {obj}')
            # sw.starlight_launcher(gridFileName, starlightFolder)
            # print('\n-Starlight finished succesfully ended')

            # Read output data
            Input_Wavelength, Input_Flux, Output_Flux, fit_output = sw.load_starlight_output(starlightFolder/'Output'/outputFile)
            z_gp = obsData['sample_data']['z_array'][i]
            Mcor, Mint = fit_output['Mcor_tot'], fit_output['Mini_tot']
            mass_galaxy = computeSSP_galaxy_mass(Mcor, 1, z_gp)
            massProcess_galaxy = computeSSP_galaxy_mass(Mint, 1, z_gp)

            if ext == ('_BR'):
                objPlotLabel = f'{obj} full spectrum'
            else:
                objPlotLabel = f'{obj} blue spectrum'


            # Plot the results
            sw.population_fraction_plots(fit_output, objPlotLabel, 'Mass_fraction', massFracPlotFile, mass_galaxy=mass_galaxy)
            sw.population_fraction_plots(fit_output, objPlotLabel, 'Light_fraction', LightFracPlotFile)
            sw.stellar_fit_comparison_plot(objPlotLabel, Input_Wavelength, Input_Flux, Output_Flux, stellarPlotFile)
            print(f'-- Mass fraction {massFracPlotFile}')

            # labelsDict = {'xlabel': r'Wavelength $(\AA)$',
            #               'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
            #               'title': f'Galaxy {objName} stellar continuum fitting'}
            #
            # idcs_plot = Input_Flux > 0.0
            #
            # fig, ax = plt.subplots(figsize=(12, 8))
            # ax.plot(wave, flux, label='--- Observed fits file')
            # ax.plot(Input_Wavelength[idcs_plot], Input_Flux[idcs_plot], label='--- Input starlight flux')
            # ax.plot(Input_Wavelength[idcs_plot], Output_Flux[idcs_plot], label='--- Output starlight fitting')
            # ax.update(labelsDict)
            # ax.legend()
            # ax.set_yscale('log')
            # plt.show()


