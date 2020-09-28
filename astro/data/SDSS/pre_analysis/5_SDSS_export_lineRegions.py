import os
import numpy as np
import pandas as pd
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from pathlib import Path

COLUMNS_TO_CLEAR = ['ion', 'intg_flux', 'intg_err', 'gauss_flux', 'gauss_err', 'eqw', 'eqw_err', 'm_continuum',
                    'n_continuum', 'std_continuum', 'amp', 'mu', 'sigma', 'amp_err', 'mu_err', 'sigma_err',
                    'pynebCode', 'pynebLabel', 'lineType','latexLabel', 'blended', 'observation', 'comments']


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
addressList = list(data_folder / file for file in fileList)

defaultExcludeLinesList = ['He1_7281A', 'He1_4388A', 'Ar4_7170A', 'H1_6563A', 'N2_6548A', 'N2_6584A']
objWithOII7300_lines = ['J003601+003307.fits', 'J012217+052044.fits', 'J012910+145935.fits']

# # Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Open lineslog
    fitsFolder, fitsFile = file_address.parent, file_address.name
    lineLogFolder, lineLogFile = fitsFolder / 'pre_analysis', fitsFile.replace('.fits', '_rawlinesLog.txt')
    masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')
    objName = fitsFile.replace('.fits', '')

    # Load line measurer object
    lm = sr.LineMesurer(linesDF_address=lineLogFolder / lineLogFile)

    # Rename regions
    lm.linesDF.rename(index={'O2_7319A_b': 'O2_7319A_m'}, inplace=True)
    lm.linesDF.rename(index={'O2_7319A': 'O2_7319A_m'}, inplace=True)
    lm.linesDF.rename(index={'He1r_4713A': 'He1_4713A'}, inplace=True)

    # Rename the labels
    for linelabel in obsData['default_line_fitting']:
        if 'A_m' in linelabel:
            line_i = linelabel[:-2]
            if line_i in lm.linesDF.index:
                lm.linesDF.rename(index={line_i: linelabel}, inplace=True)

    # Clear fit columns
    for column in COLUMNS_TO_CLEAR:
        # lm.linesDF[column] = np.nan
        lm.linesDF.drop(column, axis=1, inplace=True)

    # Save dataframe to text file
    if fitsFile in objWithOII7300_lines:
        exclude_lines = defaultExcludeLinesList
    else:
        exclude_lines = ['O2_7319A', 'O2_7330A'] + defaultExcludeLinesList
    idcs_lines = ~lm.linesDF.index.isin(exclude_lines)
    lm.save_lineslog(lm.linesDF.loc[idcs_lines], masksFolder/masksFile)