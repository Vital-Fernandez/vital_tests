import numpy as np
import pyneb as pn
import pandas as pd

from pathlib import Path
from matplotlib import rcParams, pyplot as plt, cm
from astro.papers.muse_CGCG007.muse_CGCG007_methods import read_lines_fits, compute_cHbeta, safe_cfg
from lime.tools import label_decomposition
from lime import load_cfg, save_line_log

# Cfg file
cfg_file = Path(r'D:\Pycharm Projects\vital_tests\astro\papers\muse_CGCG007\muse_CGCG007.ini')
obsCfg = load_cfg(cfg_file)

# Reddening parameters
red_curve = obsCfg['Extinction']['red_law']
R_v = obsCfg['Extinction']['R_v']

# Fits files
arms_dict = {'blue': Path(r'D:\Dropbox\Astrophysics\Data\CGCG0707_mike\AllLinesBlueArm.fits'),
             'red': Path(r'D:\Dropbox\Astrophysics\Data\CGCG0707_mike\AllLinesRedArm.fits')}

# Number of lines per fit
nights_range = range(1, 4)

# Convert fits files to dictionary of data frames
logs_dict = read_lines_fits(arms_dict, night_list=nights_range)

# Compute the reddening coefficient
for arm in list(arms_dict.keys()):
    for i_night in nights_range:

        # Get the lines for the analysis
        complete_df = logs_dict[arm][i_night]
        ion_array, wavelength_array, latexLabel_array = label_decomposition(complete_df.index.values)

        complete_df['ion'] = ion_array
        complete_df['wavelength'] = wavelength_array
        complete_df['latex_label'] = latexLabel_array

        # Find duplicated
        idcs_dup = complete_df['wavelength'].duplicated(keep=False)
        wavelengths_list = np.unique(complete_df.loc[idcs_dup].wavelength.values)

        # Add up every component
        for wave in wavelengths_list:
            idcs_comp = complete_df['wavelength'] == wave
            lineLabel = complete_df.loc[idcs_comp].index.values[0]
            flux_array, err_array = complete_df.loc[idcs_comp, 'gauss_flux'].values, complete_df.loc[idcs_comp, 'gauss_err'].values

            label_root = lineLabel[0:lineLabel.rfind('_')]
            flux, err = flux_array.sum(), np.sqrt(np.sum(np.square(err_array)))
            complete_df.loc[label_root, 'gauss_flux':'gauss_err'] = flux, err
            complete_df.loc[label_root, 'wavelength'] = wave


# Lines for the chemical analysis
lineArm_dict = {'blue': np.array(['H1_4101A', 'Ne3_3967A', 'O2_3726A', 'O2_3729A', 'O3_4363A', 'Fe3_4658A', 'He1_4471A', 'He2_4686A', 'Ar4_4711A']),
                'red':  np.array(['H1_4861A', 'O3_4959A', 'O3_5007A', 'He1_5016A', 'N2_5754A', 'He1_5876A', 'S3_6312A', 'N2_6548A', 'H1_6563A', 'N2_6584A', 'He1_6678A',
                               'S2_6716A', 'S2_6731A', 'Ar3_7135A', 'Ar3_7751A', 'H1_8750A', 'H1_8862A', 'H1_9014A', 'H1_8665A', 'H1_8598A', 'S3_9069A'])}

output_folder = Path(r'D:\Dropbox\Astrophysics\Data\CGCG0707_mike')
for i_night in nights_range:
    combined_DF = pd.DataFrame(columns=['wavelength', 'ion', 'intg_flux', 'intg_err'])
    for color, lines_list in lineArm_dict.items():
        complete_df = logs_dict[color][i_night]
        for lineLabel in lines_list:
            combined_DF.loc[lineLabel, 'wavelength'] = complete_df.loc[lineLabel, 'wavelength']
            combined_DF.loc[lineLabel, 'intg_flux'] = complete_df.loc[lineLabel, 'gauss_flux']
            combined_DF.loc[lineLabel, 'intg_err'] = complete_df.loc[lineLabel, 'gauss_err']

    # Sorte by wavelength
    combined_DF.sort_values(by=['wavelength'], inplace=True)
    ion_array, wave_array, latexLabel_array = label_decomposition(combined_DF.index.values)

    # Normalized by Hbeta
    Hbeta_flux = combined_DF.loc['H1_4861A', 'intg_flux']
    combined_DF['intg_flux'] = combined_DF['intg_flux'] / Hbeta_flux
    combined_DF['intg_err'] = combined_DF['intg_err'] / Hbeta_flux

    print(','.join(list(combined_DF.index.values)))

    # Add ions to log
    ion_array, wave_array, latexLabel_array = label_decomposition(combined_DF.index.values)
    combined_DF['ion'] = ion_array

    # Save them individually
    log_path = output_folder/f'MIKE_CGCG007_linelog_night{i_night}'
    save_line_log(combined_DF, log_path)



'''
Loop through the lineArm dict and put it there (linelabel, flux and err)

Sort it by wavelength

normalized flux by Hbeta

save it as a lines log


'''
