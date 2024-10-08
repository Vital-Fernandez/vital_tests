import numpy as np
from pathlib import Path
import src.specsiser as sr
import atpy
from src.specsiser.components.starContinuum_functions import SSPsynthesizer, computeSSP_galaxy_mass
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astro.ext_lib.starlight.plotstarlightfits import plot_fits_and_SFH
import astro.ext_lib.starlight

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

    ext = '_BR'

    # Declare files location
    fits_file = dataFolder/f'{obj}{ext}.fits'
    lineLog_file = outputFolder/f'{obj}{ext}_linesLog.txt'
    nebFluxNoNebCompFile = outputFolder/f'{obj}{ext}_obs_RemoveNebularComp.txt'
    nebCompFile = outputFolder/f'{obj}{ext}_NebFlux.txt'
    starlightOutput = starlightFolder/'Output'/f'{obj}{ext}.slOutput'
    objGasSpectrumFile = outputFolder/f'{obj}{ext}_gasSpectrum.txt'
    LightFracPlotFile = outputFolder/f'{obj}{ext}_SSP_LightFrac.png'
    stellarPlotFile = outputFolder/f'{obj}{ext}_stellarFit.png'
    specCompPlot = outputFolder/f'{obj}{ext}_specComponents.png'

    # Set wavelength and flux
    print(f'\n-- Treating {counter} :{obj}{ext}.fits')
    wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
    wave_rest = wave / (1 + z)
    idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])
    flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array

    obsWave, obsFlux = wave_rest[idx_wave], flux[idx_wave]
    nebWave, nebFlux = np.loadtxt(nebCompFile, unpack=True)
    obsNoNebWave, obsNoNebFlux = np.loadtxt(nebFluxNoNebCompFile, unpack=True)
    sw = SSPsynthesizer()
    stellarWave, inFlux, stellarFlux, fit_output = sw.load_starlight_output(starlightOutput)


    # Increase the range of Wave_S so it is greater than the observational range
    Wave_StellarExtension = np.linspace(3000.0, 3399.0, 200)
    Int_StellarExtension = np.zeros(len(Wave_StellarExtension))

    # Increase the range of Wave_S so it is greater than the observational range
    Wave_S = np.hstack((Wave_StellarExtension, stellarWave))
    Int_S = np.hstack((Int_StellarExtension, stellarFlux))

    # Resampling stellar spectra
    Interpolation = interp1d(Wave_S, Int_S, kind = 'slinear')
    Int_Stellar_Resampled = Interpolation(obsWave)

    # Plot the data
    labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                  'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                  'title': f'Galaxy {obj} HeII spectrum region'}

    idcs_plot = inFlux > 0.0
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(obsWave, obsFlux, label='Observed spectrum')
    # ax.plot(obsWave, nebFlux, label='Nebular component')
    ax.plot(obsWave, Int_Stellar_Resampled, label='Stellar component')
    ax.plot(obsWave, Int_Stellar_Resampled + nebFlux, label='Nebular + Stellar fit', linestyle=':')
    # ax.plot(obsWave, obsFluxNoStellar, label='Emission spectrum')
    # ax.plot(stellarWave[idcs_plot], inFlux[idcs_plot], label='Input starlight flux')
    # ax.plot(stellarWave[idcs_plot], stellarFlux[idcs_plot], label='Output starlight fitting')
    lineRegions = np.array([4550, 4820])
    idcs_regions = np.searchsorted(obsWave, lineRegions)
    flux_region = obsFlux[idcs_regions[0]:idcs_regions[1]]
    ax.set_xlim(lineRegions)
    ax.set_ylim(bottom=0.0, top=flux_region.mean()*1.5)

    ax.update(labelsDict)
    ax.legend()
    # ax.set_yscale('log')
    plt.show()
    # plt.savefig(lineLogFolder/specCompPlot, bbox_inches='tight')
