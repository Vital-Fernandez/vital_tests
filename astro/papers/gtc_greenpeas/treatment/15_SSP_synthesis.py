import numpy as np
import os
from pathlib import Path
import src.specsiser as sr
from src.specsiser.physical_model.starContinuum_functions import SSPsynthesizer, computeSSP_galaxy_mass
from scipy.interpolate import interp1d
from astro.papers.gtc_greenpeas.common_methods import double_arm_redCorr
import matplotlib.pyplot as plt

if os.name != 'nt':
    conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
    dataFolder = Path('/home/vital/Dropbox/Astrophysics/Papers/gtc_greenpeas/data')
    resultsFolder = Path('/home/vital/Dropbox/Astrophysics/Papers/gtc_greenpeas/treatment')
    starlight_folder = Path('/home/vital/Dropbox/Astrophysics/Tools/Starlight')
else:
    conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
    dataFolder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data')
    resultsFolder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/treatment')
    starlight_folder = Path('D:/Dropbox/Astrophysics/Tools/Starlight')

obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)
objList = obsData['file_information']['object_list']

fileList = obsData['file_information']['files_list']
idx_band = int(obsData['file_information']['band_flux'])

z_array = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
arm_wave_boundary = obsData['sample_data']['w_div']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

counter = 0
ext = '_BR'
cycle = 'c3'
cycle_ref = 'Second_cycle'

for i, obj in enumerate(objList):

    # if i == 1:

        labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                      'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                      'title': f'Galaxy {obj}{ext} nebular continuum calculation'}

        # Declare files location
        fits_file = dataFolder / f'{obj}{ext}.fits'
        objFolder = resultsFolder / f'{obj}'
        lineLog_file = objFolder / f'{obj}{ext}_linesLog.txt'
        results_file = objFolder / f'{obj}{ext}_measurements.txt'
        objMask = objFolder / f'{obj}{ext}_mask.txt'
        nebCompFile = objFolder/f'{obj}{ext}_NebFlux_{cycle}.txt'

        # Names for new data
        run_ref = f'{obj}{ext}_{cycle}'
        massFracPlotFile = objFolder / f'{obj}{ext}_SSP_MasFrac_{cycle}.png'
        LightFracPlotFile = objFolder / f'{obj}{ext}_SSP_LightFrac_{cycle}.png'
        stellarPlotFile = objFolder / f'{obj}{ext}_stellarFit_{cycle}.png'
        maskFile = starlight_folder/'Masks'/f'{obj}{ext}_{cycle}_Mask.lineslog'
        maskPlotFile = objFolder / f'{obj}{ext}_maskAndFlags_{cycle}.png'
        stellarFluxFile = objFolder / f'{obj}{ext}_stellarFlux_{cycle}.txt'

        results_dict = sr.loadConfData(results_file, group_variables=False)

        print(f'\n-- Treating: {obj}{ext}.fits')
        wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
        lm = sr.LineMesurer(wave, flux, redshift=z_array[i], crop_waves=(wmin_array[i], wmax_array[i]))

        # Recover reddening correction
        cHbeta = results_dict['Initial_values']['cHbeta_BR_Hbeta_Hgamma_Hdelta']
        int_spec, corr_spec = double_arm_redCorr(lm.wave, lm.flux, arm_wave_boundary[i], red_law, RV, cHbeta)

        # Add new masks
        linesDF = sr.lineslogFile_to_DF(lineLog_file)
        ini_mask, end_points = obsData[obj]['ini_mask_array'], obsData[obj]['end_mask_array']
        labels = ['cont_mask_' + str(int(x)) for x in ini_mask]
        for j, label_mask in enumerate(labels):
            linesDF.loc[labels[j], ['w3', 'w4']] = ini_mask[j], end_points[j]

        # Load spectra
        print(f'\n-- Treating: {obj}{ext}.fits')

        nebWave, nebFlux = np.loadtxt(nebCompFile, unpack=True)
        specFlux = lm.flux - nebFlux

        # Starlight wrapper
        sw = SSPsynthesizer()

        # Generate starlight files
        idcs_lines = ~linesDF.index.str.contains('_b')
        gridFileName, outputFile, saveFolder, waveResample, fluxResample = sw.generate_starlight_files(starlight_folder,
                                                                                                       run_ref,
                                                                                                       lm.wave,
                                                                                                       specFlux,
                                                                                                       linesDF.loc[idcs_lines])

        if os.name != 'nt':
            # Launch starlight
            print(f'\n-Initiating starlight: {obj}')
            sw.starlight_launcher(gridFileName, starlight_folder)
            print('\n-Starlight finished succesfully ended')

        # Read output data
        stellar_Wave, obj_input_flux, stellar_flux, fit_output = sw.load_starlight_output(saveFolder/outputFile)
        z_gp = obsData['sample_data']['z_array'][i]
        Mcor, Mint = fit_output['Mcor_tot'], fit_output['Mini_tot']
        mass_galaxy = computeSSP_galaxy_mass(Mcor, 1, z_gp)
        massProcess_galaxy = computeSSP_galaxy_mass(Mint, 1, z_gp)

        # Plot the results
        plot_label = f'{obj} spectrum' if ext == '_BR' else f'{obj} blue arm spectrum'
        sw.population_fraction_plots(fit_output, plot_label, 'Mass_fraction', massFracPlotFile, mass_galaxy=mass_galaxy)
        sw.population_fraction_plots(fit_output, plot_label, 'Light_fraction', LightFracPlotFile)
        sw.stellar_fit_comparison_plot(plot_label, stellar_Wave, obj_input_flux, stellar_flux, stellarPlotFile)
        sw.mask_plot(fit_output, obj, lm.wave, specFlux, stellar_Wave, obj_input_flux, maskFile, maskPlotFile)

        # Store starlight configuration values for linux run
        starlight_cfg = {'gridFileName': gridFileName,
                         'outputFile': outputFile,
                         'saveFolder': saveFolder.as_posix(),
                         'Galaxy_mass': mass_galaxy,
                         'Chi2': fit_output['Chi2']}
        sr.parseConfDict(results_file, starlight_cfg, f'Starlight_run_{cycle}', clear_section=True)

        # ----- Save the object flus without the stellar component
        #Increase the range of Wave_S so it is greater than the observational range
        Wave_StellarExtension = np.linspace(3000.0, 3399.0, 200)
        Int_StellarExtension = np.zeros(len(Wave_StellarExtension))

        #Increase the range of Wave_S so it is greater than the observational range
        Wave_S = np.hstack((Wave_StellarExtension, stellar_Wave))
        Int_S = np.hstack((Int_StellarExtension, stellar_flux))

        #Resampling stellar spectra
        Interpolation = interp1d(Wave_S, Int_S, kind = 'slinear')
        flux_Stellar_Resampled = Interpolation(lm.wave)

        # Save the non object spectrum without stellar component
        stellarFlux = flux_Stellar_Resampled
        np.savetxt(stellarFluxFile, np.transpose(np.array([lm.wave, stellarFlux])), fmt="%7.1f %10.4e")

        # compare adding componennts
        wave_neb, flux_neb = np.loadtxt(nebCompFile, unpack=True)
        wave_star, flux_star = np.loadtxt(stellarFluxFile, unpack=True)

        # # Plot spectra components
        # fig, ax = plt.subplots(figsize=(12, 8))
        # ax.plot(lm.wave, lm.flux, label='Object flux')
        # # ax.plot(stellar_Wave, obj_input_flux, label='Input starlight spectrum', linestyle='--')
        # ax.plot(wave_neb, flux_neb, label='Nebular flux')
        # ax.plot(wave_star, flux_star, label='Stellar flux')
        # ax.plot(wave_star, flux_star + flux_neb, label='Combined continuum', linestyle=':')
        # ax.legend()
        # ax.set_yscale('log')
        # plt.tight_layout()
        # plt.show()
