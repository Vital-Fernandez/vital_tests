from pathlib import Path
import numpy as np
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.starContinuum_functions import StarlightWrapper
import matplotlib.pyplot as plt

# Import the observation data
obsData = sr.loadConfData('../gtc_greenpeas_data.ini', group_variables=False)
# starlightFolder = Path('/home/vital/Astro-data/osiris-Ricardo/starlight')
starlightFolder = Path('D:/Google drive/Astrophysics/Datos/osiris-Ricardo/starlight')
data_folder = Path(obsData['file_information']['data_folder'])
file_list = obsData['file_information']['files_list']
addressList = list(data_folder/file for file in file_list)
flux_norm = obsData['sample_data']['norm_flux']

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Establish files location
    objName = obsData['file_information']['object_list'][i]
    fitsFolder, fitsFile = file_address.parent, file_address.name
    masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')
    lineLogFolder = fitsFolder/'flux_analysis'
    lineLogFile, stellarPlotFile = fitsFile.replace('.fits', '_linesLog.txt'), fitsFile.replace('.fits', '_stellarFit.png')
    nebFluxNoNebCompFile = fitsFile.replace('.fits', '_obs_RemoveNebularComp.txt')
    print(f'\n-{objName}')

    # Get fits data
    wave, flux, header = sr.import_fits_data(file_address, instrument='OSIRIS')
    z_mean = obsData['sample_data']['z_array'][i]
    wmin_array, wmax_array = obsData['sample_data']['wmin_array'], obsData['sample_data']['wmax_array']

    # Define wave and flux ranges
    wave_rest = wave / (1 + z_mean)
    idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])

    # Load the data
    obsWave, obsFlux = wave_rest[idx_wave], flux[idx_wave]
    specWave, specFlux = np.loadtxt(lineLogFolder/nebFluxNoNebCompFile, unpack=True)

    # Measuring objects
    lm = sr.LineMesurerGUI(specWave, specFlux, lineLogFolder/lineLogFile, normFlux=flux_norm)
    sw = StarlightWrapper()

    # Generate starlight files
    gridFileName, outputFile, outputFolder, waveResample, fluxResample = sw.generate_starlight_files(starlightFolder,
                                                                                                     objName,
                                                                                                     specWave,
                                                                                                     specFlux,
                                                                                                     lm.linesDF)

        # Launch starlight
        # print(f'\n-Initiating starlight: {objName}')
        # sw.starlight_launcher(gridFileName, starlightFolder)
        # print('\n-Starlight finished succesfully ended')

        # Read output data
        # Input_Wavelength, Input_Flux, Output_Flux, fit_output = sw.load_starlight_output(outputFolder/outputFile)

        # Plot the results
        # sw.population_fraction_plots(fit_output, objName, 'Mass_fraction', lineLogFolder/f'{objName}_SSP_MasFrac.png')
        # sw.population_fraction_plots(fit_output, objName, 'Light_fraction', lineLogFolder/f'{objName}_SSP_LightFrac.png')
        # sw.stellar_fit_comparison_plot(objName, Input_Wavelength, Input_Flux, Output_Flux, lineLogFolder/stellarPlotFile)



