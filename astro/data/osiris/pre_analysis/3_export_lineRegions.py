import lmfit
import numpy as np
import pandas as pd
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from pathlib import Path

COLUMNS_TO_CLEAR = ['intg_flux' ,'intg_err','gauss_flux','gauss_err','eqw','eqw_err', 'm_continuum', 'n_continuum',
                   'std_continuum','amp','mu','sigma','amp_err','mu_err','sigma_err']

# Import the observation data
obsData = sr.loadConfData('../gtc_greenpeas_data.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
data_folder = Path(obsData['file_information']['data_folder'])
file_list = obsData['file_information']['files_list']
addressList = list(data_folder/file for file in file_list)

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Establish files location
    objName = obsData['file_information']['object_list'][i]
    fitsFolder, fitsFile = file_address.parent, file_address.name
    logFolder, logFile = fitsFolder/'pre_analysis', fitsFile.replace('.fits', '_rawlinesLog.txt')
    masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')

    # Load line measurer object
    lm = sr.LineMesurer(linesDF_address=logFolder / logFile)
    print(lm.label_formatter('O3_5007A'))

    # Special corrections
    if 'He1r_4713A' in lm.linesDF.index:
        lm.linesDF.rename(index={'He1r_4713A': 'Ar4_4711A'}, inplace=True)
        lm.linesDF.loc['Ar4_4711A'] = sr._linesDb.loc['Ar4_4711A']
    if 'H1_3704A' in lm.linesDF.index:
        lm.linesDF.rename(index={'H1_3704A': 'He1_3704A'}, inplace=True)
        lm.linesDF.loc['He1_3704A'] = sr._linesDb.loc['He1_3704A']
    if 'H1_3704A_b' in lm.linesDF.index:
        lm.linesDF.rename(index={'H1_3704A_b': 'He1_3704A'}, inplace=True)
        lm.linesDF.loc['He1_3704A'] = sr._linesDb.loc['He1_3704A']
    if 'Ar3_7751A' in lm.linesDF.index:
        lm.linesDF.loc['Ar3_7751A', 'latexLabel'] = '$7751\,[ArIII]$'

    # Label the blended lines
    default_blended_groups = obsData['blended_groups']
    blended_list = list(default_blended_groups.keys())
    conversion_dict = dict(zip(map(lambda each: each.strip('_b'), blended_list), blended_list))
    lm.linesDF.rename(index=conversion_dict, inplace=True)

    for line in default_blended_groups:
        if line in obsData[f'{objName}_blended_lines']:
            line_components = obsData[f'{objName}_blended_lines'][line]
        else:
            line_components = default_blended_groups[line]
        lm.linesDF.loc[line, 'blended'] = line_components
        for component in line_components.split('-'):
            if component in lm.linesDF.index:
                lm.linesDF.drop(component, inplace=True)
        # print('blended', line, default_blended_groups[line], lm.label_formatter(default_blended_groups[line]))
        lm.linesDF.loc[line, 'latexLabel'] = lm.label_formatter(default_blended_groups[line])

    # Label the merged lines
    default_merged_groups = obsData['merged_groups']
    merged_list = list(default_merged_groups.keys())
    conversion_dict = dict(zip(map(lambda each: each.strip('_m'), merged_list), merged_list))
    lm.linesDF.rename(index=conversion_dict, inplace=True)

    for line in default_merged_groups:
        if line in lm.linesDF.index:
            lm.linesDF.loc[line, 'blended'] = default_merged_groups[line]
        components_list = default_merged_groups[line]
        for component in default_merged_groups:
            if component in lm.linesDF.index:
                lm.linesDF.drop(component)
        lm.linesDF.loc[line, 'latexLabel'] = lm.label_formatter(default_merged_groups[line])

    # Clear fit columns
    for column in COLUMNS_TO_CLEAR:
        lm.linesDF[column] = np.nan

    # Set all lines to true detection
    lm.linesDF['observation'] = 'detected'

    # New column continuum
    lm.linesDF.insert(loc=19, column='cont', value=np.nan)

    # Save dataframe to text file
    lm.save_lineslog(lm.linesDF, masksFolder/masksFile)

