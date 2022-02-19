import numpy as np
import os
from pathlib import Path
import src.specsiser as sr
from src.specsiser.components.starContinuum_functions import SSPsynthesizer, computeSSP_galaxy_mass
from scipy.interpolate import interp1d
from astro.papers.gtc_greenpeas.common_methods import double_arm_redCorr
import pyneb as pn

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

w_div_array = obsData['sample_data']['w_div']
red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

ext = 'BR'
cycle = 'it1'

for i, obj in enumerate(objList):

    print(f'\n-- Treating: {obj}{ext}.fits')

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    fits_file = dataFolder / f'{obj}_{ext}.fits'
    objMask = dataFolder / f'{obj}_{ext}_mask.txt'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
    nebCompFile = objFolder/f'{obj}_{ext}_nebFlux_{cycle}.txt'

    # Declare output files
    massFracPlotFile = objFolder / f'{obj}_{ext}_SSP_MasFrac_{cycle}.png'
    LightFracPlotFile = objFolder / f'{obj}_{ext}_SSP_LightFrac_{cycle}.png'
    stellarPlotFile = objFolder / f'{obj}_{ext}_stellarFit_{cycle}.png'
    maskFile = starlight_folder/'Masks'/f'{obj}_{ext}_{cycle}_Mask.lineslog'
    maskPlotFile = objFolder / f'{obj}_{ext}_maskAndFlags_{cycle}.png'
    stellarFluxFile = objFolder / f'{obj}_{ext}_stellarFlux_{cycle}.txt'

    # Load data
    results_dict = sr.loadConfData(results_file, group_variables=False)

    linesDF = sr.lineslogFile_to_DF(lineLog_file)

    wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
    flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
    lm = sr.LineMesurer(wave, flux, redshift=z_array[i], crop_waves=(wmin_array[i], wmax_array[i]))

    nebWave, nebFlux = np.loadtxt(nebCompFile, unpack=True)

    # Spectrum extinction correction
    cHbeta_label = obsData[obj]['cHbeta_label']
    cHbeta = np.array(results_dict[f'Extinction_{cycle}'][cHbeta_label], dtype=float)
    int_spec, corr_spec = double_arm_redCorr(lm.wave, lm.flux, w_div_array[i], red_law, RV, cHbeta)

    # Add new masks
    linesDF = sr.lineslogFile_to_DF(lineLog_file)
    ini_mask, end_points = obsData[obj]['ini_mask_array'], obsData[obj]['end_mask_array']
    labels = ['cont_mask_' + str(int(x)) for x in ini_mask]
    for j, label_mask in enumerate(labels):
        linesDF.loc[labels[j], ['w3', 'w4']] = ini_mask[j], end_points[j]

    # Remove nebular component spectra
    specFlux = lm.flux - nebFlux

    # Starlight wrapper
    sw = SSPsynthesizer()

    # Generate starlight files
    clip_value = obsData[obj]['starlight_clipping']
    idcs_lines = ~linesDF.index.str.contains('_b')
    gridFileName, outputFile, saveFolder, waveResample, fluxResample = sw.generate_starlight_files(starlight_folder,
                                                                                                   f'{obj}_{ext}_{cycle}',
                                                                                                   lm.wave,
                                                                                                   specFlux,
                                                                                                   linesDF.loc[idcs_lines],
                                                                                                   clip_value=clip_value)

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
    idcs_below_20Myr = fit_output['DF'].age_j < 2*10**7
    mass_galaxy_20Myr_percent = np.sum(fit_output['DF'].loc[idcs_below_20Myr, 'Mcor_j'].values)

    # Store starlight configuration values for linux run
    rc = pn.RedCorr(R_V=RV, E_BV=fit_output['Av_min']/RV, law=red_law)
    cHbeta_star = rc.cHbetaFromEbv(fit_output['Av_min']/RV)
    starlight_cfg = {'gridFileName': gridFileName,
                     'outputFile': outputFile,
                     'saveFolder': saveFolder.as_posix(),
                     'Galaxy_mass_Current': mass_galaxy,
                     'Galaxy_mass_Prosessed': massProcess_galaxy,
                     'Galaxy_mass_Percentbelow20Myr': mass_galaxy_20Myr_percent,
                     'Chi2': fit_output['Chi2'],
                     'A_V_stellarr': fit_output['Av_min'],
                     'cHbeta_stellar': cHbeta_star,
                     'PixelMeanDevPer': fit_output['SumXdev'],
                     'SN': fit_output['SignalToNoise_magnitudeWave']}
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

    # Plot the results
    plot_label = f'{obj} spectrum' if ext == '_BR' else f'{obj} blue arm spectrum'
    sw.population_fraction_plots(fit_output, plot_label, 'Mass_fraction', massFracPlotFile, mass_galaxy=mass_galaxy)
    sw.population_fraction_plots(fit_output, plot_label, 'Light_fraction', LightFracPlotFile)
    sw.stellar_fit_comparison_plot(f'{obj}_{ext}_{cycle}', lm.wave, lm.flux, nebCompFile, stellarFluxFile, stellarPlotFile)
    sw.mask_plot(fit_output, obj, lm.wave, specFlux, stellar_Wave, stellar_flux, obj_input_flux, maskFile, maskPlotFile)




