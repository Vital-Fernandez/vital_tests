import os
import numpy as np
import pandas as pd
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from pathlib import Path


def list_files(directory, extension):
    output_list = []
    for file in os.listdir(directory):
        if file.endswith(extension):
            output_list.append(os.path.join(directory, file))
    return output_list


# Import the observation data
obsData = sr.loadConfData('../flux_comparison.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['file_information']['data_folder'])
fileList = list_files(data_folder, '.fits')
addressList = list(data_folder/file for file in fileList)

for file_address in addressList:

    print(file_address)

    # Set and crop the wavelength
    wave_rest, flux, header = sr.import_fits_data(file_address, instrument='SDSS')
    idx_wave = (wave_rest >= obsData['sample_data']['wmin_array']) & (wave_rest <= obsData['sample_data']['wmax_array'])

    # Analyse the spectrum
    lm = sr.LineMeasurer(wave_rest[idx_wave], flux[idx_wave])
    #lm.plot_spectrum_components()

    # Normalize
    noise_region = obsData['sample_data']['noiseRegion_array']
    norm_flux = lm.continuum_remover(noise_region)

    # Detect the observed lines
    blendedLinesList = obsData['sample_data']['blendedLines_list']
    obsLinesTable = lm.line_finder(norm_flux, noiseWaveLim=noise_region, intLineThreshold=3)
    obsLinesDF = lm.match_lines(obsLinesTable, linesDb, blendedLineList=blendedLinesList)
    # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=obsLinesDF)

    # Measure line fluxes
    idcsObsLines = (obsLinesDF.observation == 'detected')
    obsLines = obsLinesDF.loc[idcsObsLines].index.values
    blended_groups = obsData['blended_groups']

    for lineLabel in obsLines:
        print(f'- {lineLabel}:')

        # Declare regions data
        wave_regions = obsLinesDF.loc[lineLabel, 'w1':'w6'].values
        idcsLinePeak, idcsContinua = lm.define_masks(wave_regions)

        # Identify line regions
        lm.line_properties(idcsLinePeak, idcsContinua, bootstrap_size=1000)

        # Perform gaussian fitting
        lm.linesDF = obsLinesDF
        lm.line_fit('lmfit', lineLabel, idcsLinePeak, idcsContinua, blended_groups=blended_groups)

        # Store results in database
        lm.results_to_database(lineLabel, obsLinesDF)

    # Save dataframe to text file
    linesLogAddress = str(file_address).replace('.fits', '_rawlinesLog.txt')
    lm.save_lineslog(obsLinesDF.loc[idcsObsLines], linesLogAddress)
