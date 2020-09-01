from pathlib import Path
import numpy as np
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.starContinuum_functions import StarlightWrapper
import matplotlib.pyplot as plt

# Import the observation data
obsData = sr.loadConfData('../gtc_greenpeas_data.ini', group_variables=False)
starlightFolder = Path('/home/vital/Astro-data/osiris-Ricardo/starlight')
# starlightFolder = Path('D:/Google drive/Astrophysics/Datos/osiris-Ricardo/starlight')
data_folder = Path(obsData['file_information']['data_folder'])
file_list = obsData['file_information']['files_list']
addressList = list(data_folder/file for file in file_list)
flux_norm = obsData['sample_data']['norm_flux']

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    if i == 2:

        # Establish files location
        objName = obsData['file_information']['object_list'][i]
        fitsFolder, fitsFile = file_address.parent, file_address.name
        masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')
        lineLogFolder = fitsFolder/'flux_analysis'
        lineLogFile, stellarPlotFile = fitsFile.replace('.fits', '_linesLog.txt'), fitsFile.replace('.fits', '_stellarComp.png')
        nebFluxNoNebCompFile = fitsFile.replace('.fits', '_NoNebFlux.txt')
        print(f'-{objName}')

        # Load the data
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

        # # Launch starlight
        # print(f'\n-Initiating starlight: {objName}')
        # sw.starlight_launcher(gridFileName, starlightFolder)
        # print('\n-Starlight finished succesfully ended')

        # Read output data
        Input_Wavelength, Input_Flux, Output_Flux, Parameters = sw.load_starlight_output(outputFolder/outputFile)

        maskPix, clipPix, flagPix = Parameters['MaskPixels'], Parameters['ClippedPixels'], Parameters['FlagPixels']

        # Plot spectra components
        labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                      'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                      'title': f'Galaxy {objName} stellar continuum fitting'}

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.plot(Input_Wavelength, Input_Flux, label='Input starlight flux')
        ax.plot(Input_Wavelength, Output_Flux, label='Output starlight fitting')
        ax.update(labelsDict)
        ax.legend()
        ax.set_yscale('log')
        # plt.savefig(lineLogFolder/nebPlotFile, bbox_inches='tight')
        plt.show()
        fig.clear()


