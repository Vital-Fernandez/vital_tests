from pathlib import Path
import numpy as np
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.starContinuum_functions import StarlightWrapper
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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

    sw = StarlightWrapper()

    # Establish files location
    objName = obsData['file_information']['object_list'][i]
    fitsFolder, fitsFile = file_address.parent, file_address.name
    masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')
    lineLogFolder = fitsFolder/'flux_analysis'
    lineLogFile = fitsFile.replace('.fits', '_linesLog.txt')
    specCompPlot = fitsFile.replace('.fits', '_specComponents.png')
    nebFluxNoNebCompFile = fitsFile.replace('.fits', '_obs_RemoveNebularComp.txt')
    nebCompFile = fitsFile.replace('.fits', '_NebFlux.txt')
    outputFileAddress = starlightFolder/f'Output/{objName}.slOutput'
    objGasSpectrumFile = fitsFolder/'flux_analysis'/f'{objName}_gasSpectrum.txt'
    print(f'-{objName}')

    # Get fits data
    wave, flux, header = sr.import_fits_data(file_address, instrument='OSIRIS')
    z_mean = obsData['sample_data']['z_array'][i]
    wmin_array, wmax_array = obsData['sample_data']['wmin_array'], obsData['sample_data']['wmax_array']

    # Define wave and flux ranges
    wave_rest = wave / (1 + z_mean)
    idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])

    # Load the data
    obsWave, obsFlux = wave_rest[idx_wave], flux[idx_wave]
    nebWave, nebFlux = np.loadtxt(lineLogFolder/nebCompFile, unpack=True)
    obsNoNebWave, obsNoNebFlux = np.loadtxt(lineLogFolder/nebFluxNoNebCompFile, unpack=True)
    stellarWave, inFlux, stellarFlux, fit_output = sw.load_starlight_output(outputFileAddress)

    #Increase the range of Wave_S so it is greater than the observational range
    Wave_StellarExtension = np.linspace(3000.0,3399.0,200)
    Int_StellarExtension = np.zeros(len(Wave_StellarExtension))

    #Increase the range of Wave_S so it is greater than the observational range
    Wave_S = np.hstack((Wave_StellarExtension, stellarWave))
    Int_S = np.hstack((Int_StellarExtension, stellarFlux))

    #Resampling stellar spectra
    Interpolation = interp1d(Wave_S, Int_S, kind = 'slinear')
    Int_Stellar_Resampled = Interpolation(obsWave)

    # Save the non object spectrum without stellar component
    np.savetxt(objGasSpectrumFile, np.transpose(np.array([obsWave, obsFlux-Int_Stellar_Resampled])), fmt="%7.1f %10.4e")

    # Plot the data
    labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                  'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                  'title': f'Galaxy {objName} spectrum components comparison'}

    idcs_plot = inFlux > 0.0
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(obsWave, obsFlux, label='Observed spectrum')
    ax.plot(obsWave, nebFlux, label='Nebular component')
    ax.plot(obsWave, Int_Stellar_Resampled, label='Stellar component')
    ax.plot(obsWave, Int_Stellar_Resampled + nebFlux, label='Nebular + Stellar', linestyle=':')
    # ax.plot(obsWave, obsFlux - Int_Stellar_Resampled, label='Emission spectrum')
    # ax.plot(stellarWave[idcs_plot], inFlux[idcs_plot], label='Input starlight flux')
    # ax.plot(stellarWave[idcs_plot], stellarFlux[idcs_plot], label='Output starlight fitting')
    ax.update(labelsDict)
    ax.legend()
    ax.set_yscale('log')
    plt.savefig(lineLogFolder/specCompPlot, bbox_inches='tight')
    # plt.show()
    # fig.clear()