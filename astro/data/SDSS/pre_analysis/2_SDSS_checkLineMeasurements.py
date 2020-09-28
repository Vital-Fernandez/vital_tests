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
obsData = sr.loadConfData('D:/Pycharm Projects/vital_tests/astro/data/SDSS/flux_comparison.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['file_information']['data_folder'])
fileList = list_files(data_folder, '.fits')
addressList = list(data_folder/file for file in fileList)


# # Analyse the spectrum
for i, file_address in enumerate(addressList):

    # if 'J1026' in str(file_address):


    # Open lineslog
    fitsFolder, fitsFile = file_address.parent, file_address.name
    lineLogFolder, lineLogFile = fitsFolder / 'pre_analysis', fitsFile.replace('.fits', '_rawlinesLog.txt')

    print(fitsFile)


    # Set and crop the wavelength
    wave_rest, flux, header = sr.import_fits_data(fitsFolder/fitsFile, instrument='SDSS')
    idx_wave = (wave_rest >= obsData['sample_data']['wmin_array']) & (wave_rest <= obsData['sample_data']['wmax_array'])

    # Load line measurer object
    lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave], lineLogFolder / lineLogFile)

    # Plot the matched lines:
    lm.plot_detected_lines(lm.linesDF, ncols=5)

    # # Save dataframe to text file
    # linesLogAddressTreated = str(file_address).replace('.fits', '_treatedlinesLog.txt')
    # lm.save_lineslog(lm.linesDF)