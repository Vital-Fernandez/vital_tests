import lmfit
import numpy as np
import pandas as pd
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from pathlib import Path

# Import the observation data
obsData = sr.loadConfData('gtc_greenpeas_data.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['file_information']['data_folder'])
file_list = obsData['file_information']['files_list']
addressList = list(data_folder/file for file in file_list)

for i, file_address in enumerate(addressList):

    # Get fits data
    wave, flux, header = sr.import_fits_data(file_address, instrument='OSIRIS')

    # Get observation data
    objName = header['OBJECT']
    objReference = obsData['file_information']['object_list'][i]
    objWaves = obsData['sample_data'][f'{objReference}_obsWaves_array']

    # Compute the redshifts
    # redshifts = (objWaves/refLines) - 1
    # z_mean, z_std = redshifts.mean(), redshifts.std()
    # wave_rest = wave / (1 + z_mean)
    # print(objReference, z_mean, z_std)

    # Set and crop the wavelength
    z_mean = obsData['sample_data']['z_array'][i]
    wave_rest = wave / (1 + z_mean)
    wmin_array, wmax_array = obsData['sample_data']['wmin_array'], obsData['sample_data']['wmax_array']
    idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])

    # Analyse the spectrum
    lm = sr.LineMeasurer(wave_rest[idx_wave], flux[idx_wave])

    # Normalize
    noise_region = obsData['sample_data']['noiseRegion_array']
    norm_flux = lm.continuum_remover(noise_region)

    # Detect the observed lines
    blendedLinesList = obsData['sample_data']['blendedLines_list']
    obsLinesTable = lm.line_finder(norm_flux, noiseWaveLim=noise_region, intLineThreshold=3)
    obsLinesDF = lm.match_lines(obsLinesTable, linesDb, blendedLineList=blendedLinesList)
    # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=obsLinesDF)

    # Plot the matched lines:
    lm.plot_detected_lines(obsLinesDF)

    # Measure line fluxes
    idcsObsLines = (linesDb.observation == 'detected')
    obsLines = linesDb.loc[idcsObsLines].index.values

    for lineLabel in obsLines:
        print(f'- {lineLabel}:')

        # Declare regions data
        wave_regions = linesDb.loc[lineLabel, 'w1':'w6'].values
        idcsLinePeak, idcsContinua = lm.define_masks(wave_regions)

        # Identify line regions
        lm.line_properties(idcsLinePeak, idcsContinua, bootstrap_size=500)

        # Perform gaussian fitting
        lm.line_fitting(idcsLinePeak, idcsContinua, bootstrap_size=500)

        # Store results in database
        lm.results_to_database(lineLabel, linesDb)

    # Save dataframe to text file
    linesLogAddress = str(file_address).replace('.fits', '_rawlinesLog.txt')
    lm.save_lineslog(linesDb.loc[idcsObsLines], linesLogAddress)

    # Plot the matched lines:
    lm.plot_detected_lines(linesDb)

