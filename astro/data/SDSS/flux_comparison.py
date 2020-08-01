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

obsData = sr.loadConfData('flux_comparison.ini')
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['data_folder'])
fileList = list_files(data_folder, '.fits')
print(fileList)
addressList = list(data_folder/file for file in fileList)

for file_address in addressList:

    wave_rest, flux, header = sr.import_fits_data(file_address, instrument='SDSS')
    idx_wave = (wave_rest >= obsData['wmin_array']) & (wave_rest <= obsData['wmax_array'])

    # Analyse the spectrum
    lm = sr.LineMeasurer(wave_rest[idx_wave], flux[idx_wave])
    #lm.plot_spectrum_components()

    # Normalize
    norm_flux = lm.continuum_remover(obsData['noiseRegion_array'])

    # Detect the observed lines
    lm.blendedLinesList = obsData['noiseRegion_array']
    obsLinesTable = lm.line_finder(norm_flux, noiseWaveLim=obsData['noiseRegion_array'], intLineThreshold=3)
    obsLinesDF = lm.match_lines(obsLinesTable, linesDb, blendedLineList=obsData['blendedLines_list'])
    lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=obsLinesDF)

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
    linesLogAddress = str(file_address).replace('.fits', '_linesLog.txt')
    lm.save_lineslog(linesDb.loc[idcsObsLines], linesLogAddress)

    # Plot the matched lines:
    lm.plot_detected_lines(linesDb)
