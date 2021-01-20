import numpy as np
from pathlib import Path
import src.specsiser as sr
from src.specsiser.physical_model.starContinuum_functions import SSPsynthesizer, computeSSP_galaxy_mass
import matplotlib.pyplot as plt

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
dataFolder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data')
resultsFolder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/treatment')
starlight_folder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/starlight')

obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)
objList = obsData['file_information']['object_list']
# dataFolder = Path(obsData['file_information']['data_folder'])
# resultsFolder = Path(obsData['file_information']['results_folder'])
# starlightFolder = Path(obsData['SSP_synthesis']['starlight_folder'])

fileList = obsData['file_information']['files_list']
idx_band = int(obsData['file_information']['band_flux'])

z_array = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

counter = 0
ext = '_BR'

for i, obj in enumerate(objList):

    # Declare files location
    fits_file = dataFolder / f'{obj}{ext}.fits'
    objFolder = resultsFolder / f'{obj}'
    lineLog_file = objFolder / f'{obj}{ext}_linesLog.txt'
    results_file = objFolder / f'{obj}{ext}_measurements.txt'
    objMask = objFolder / f'{obj}{ext}_mask.txt'
    nebFluxNoNebCompFile = objFolder / f'{obj}{ext}_obs_RemoveNebularComp.txt'

    massFracPlotFile = objFolder / f'{obj}{ext}_SSP_MasFrac.png'
    LightFracPlotFile = objFolder / f'{obj}{ext}_SSP_LightFrac.png'
    stellarPlotFile = objFolder / f'{obj}{ext}_stellarFit.png'

    results_dict = sr.loadConfData(results_file, group_variables=False)
    linesDF = sr.lineslogFile_to_DF(lineLog_file)

    # Physical parameters
    cHbeta = results_dict['Initial_values']['cHbeta_BR_Hbeta_Hgamma_Hdelta']

    # Load spectra
    print(f'\n-- Treating: {obj}{ext}.fits')
    specWave, specFlux = np.loadtxt(nebFluxNoNebCompFile, unpack=True)

    # Starlight wrapper
    sw = SSPsynthesizer()

    # Generate starlight files
    run_ref = f'{obj}{ext}'
    idcs_lines = ~linesDF.index.str.contains('_b')
    gridFileName, outputFile, saveFolder, waveResample, fluxResample = sw.generate_starlight_files(starlight_folder,
                                                                                                   run_ref,
                                                                                                   specWave,
                                                                                                   specFlux,
                                                                                                   linesDF.loc[idcs_lines])

    # Store starlight configuration values for linux run
    starlight_cfg = {'gridFileName': gridFileName, 'outputFile': outputFile, 'saveFolder': saveFolder.as_posix()}
    sr.parseConfDict(results_file, starlight_cfg, 'Starlight_run1')

    # # Launch starlight
    # print(f'\n-Initiating starlight: {obj}')
    # sw.starlight_launcher(gridFileName, starlightFolder)
    # print('\n-Starlight finished succesfully ended')

    # Read output data
    Input_Wavelength, Input_Flux, Output_Flux, fit_output = sw.load_starlight_output(saveFolder/outputFile)
    z_gp = obsData['sample_data']['z_array'][i]
    Mcor, Mint = fit_output['Mcor_tot'], fit_output['Mini_tot']
    mass_galaxy = computeSSP_galaxy_mass(Mcor, 1, z_gp)
    massProcess_galaxy = computeSSP_galaxy_mass(Mint, 1, z_gp)

    # Plot the results
    plot_label = f'{obj} spectrum' if ext == '_BR' else f'{obj} blue arm spectrum'
    sw.population_fraction_plots(fit_output, plot_label, 'Mass_fraction', massFracPlotFile, mass_galaxy=mass_galaxy)
    sw.population_fraction_plots(fit_output, plot_label, 'Light_fraction', LightFracPlotFile)
    sw.stellar_fit_comparison_plot(plot_label, Input_Wavelength, Input_Flux, Output_Flux, stellarPlotFile)

