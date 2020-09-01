from pathlib import Path
import numpy as np
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.starContinuum_functions import StarlightWrapper
import matplotlib.pyplot as plt

# Import the observation data
obsData = sr.loadConfData('../gtc_greenpeas_data.ini', group_variables=False)
starlightFolder = Path('D:/Google drive/Astrophysics/Datos/osiris-Ricardo/starlight')
data_folder = Path(obsData['file_information']['data_folder'])
addressList = list(data_folder/file for file in file_list)
file_list = obsData['file_information']['files_list']
flux_norm = obsData['sample_data']['norm_flux']

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    if i == 2:

        # Establish files location
        objName = obsData['file_information']['object_list'][i]
        fitsFolder, fitsFile = file_address.parent, file_address.name
        masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')
        lineLogFolder = fitsFolder/'flux_analysis'
        lineLogFile, nebPlotFile = fitsFile.replace('.fits', '_linesLog.txt'), fitsFile.replace('.fits', '_nebComp.png')
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
        # print(f'--Initiating starlight: {objName}')
        # sw.starlight_launcher(gridFileName, starlightFolder)
        # print('--Starlight finished succesfully ended')
        #
        # # Read output data
        # Input_Wavelength, Input_Flux, Output_Flux, Parameters = sw.load_starlight_output(outputFolder/outputFile)




