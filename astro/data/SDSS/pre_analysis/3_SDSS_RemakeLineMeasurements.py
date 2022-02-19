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
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/database/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['file_information']['data_folder'])
fileList = list_files(data_folder, '.fits')
addressList = list(data_folder/file for file in fileList)
fluxNorm = obsData['sample_data']['norm_flux']

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Open lineslog
    linesLogAddress = str(file_address).replace('.fits', '_treatedlinesLog.txt')
    print(linesLogAddress)

    # Set and crop the wavelength
    wave_rest, flux, header = sr.import_fits_data(file_address, instrument='SDSS')
    idx_wave = (wave_rest >= obsData['sample_data']['wmin_array']) & (wave_rest <= obsData['sample_data']['wmax_array'])

    # Load line measurer object
    lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave] / fluxNorm, linesLogAddress)

    # Measure line fluxes
    idcs_lines = (lm.linesDF['blended_label'] == 'None') | (lm.linesDF.index.str.contains('_b'))
    obsLines = lm.linesDF.loc[idcs_lines].index.values
    blended_groups = obsData['blended_groups']

    # Loop throug the liens
    for j, lineLabel in enumerate(obsLines):

        print(f'- {lineLabel}:')

        # Declare regions data
        wave_regions = lm.linesDF.loc[lineLabel, 'w1':'w6'].values
        idcsLinePeak, idcsContinua = lm.define_masks(wave_regions)

        # Identify line regions
        lm.line_properties(idcsLinePeak, idcsContinua, bootstrap_size=1000)

        # Establish line and object fit configuration # TODO put all this tricks in an external method
        fit_conf = {}
        if 'default_blended_lines' in obsData:
            fit_conf.update(obsData['default_blended_lines'])
        if f'{lineLabel}_blended_lines' in obsData:
            fit_conf.update(obsData[f'{lineLabel}_blended_lines'])

        # Perform gaussian fitting
        lm.line_fit('lmfit', lineLabel, idcsLinePeak, idcsContinua, continuum_check=True,
                    blended_groups=blended_groups, user_conf=fit_conf)

    # Save dataframe to text file
    linesLogAddressTreated = str(file_address).replace('.fits', '_treatedlinesLog.txt')
    lm.save_lineslog(lm.linesDF, linesLogAddressTreated)

    # # Plot the matched lines:
    # lm.plot_detected_lines(lm.linesDF, ncols=5)

