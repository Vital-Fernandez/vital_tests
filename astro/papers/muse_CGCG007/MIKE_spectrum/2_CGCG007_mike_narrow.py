import numpy as np
import pyneb as pn
import pandas as pd

from astropy.io import fits
from pathlib import Path
from matplotlib import rcParams, pyplot as plt, cm
from astro.papers.muse_CGCG007.muse_CGCG007_methods import read_lines_fits, compute_cHbeta, safe_cfg
from lime import load_cfg

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

# Lines for the extinction calculation
ext_lines = {'blue': ['H1_4861A_N1', 'H1_4340A_N1', 'H1_4101A_N1', 'H1_3734A', 'H1_3750A', 'H1_3770A', 'H1_3797A', 'H1_3703A'],
             'red':  ['H1_4861A_N1', 'H1_6563A_N1', 'H1_8545A', 'H1_8598A', 'H1_8665A', 'H1_8750A', 'H1_8862A', 'H1_9014A']}

# Compute the reddening coefficient
redcorr_dict = {arm: {i_night: None for i_night in nights_range} for arm in arms_dict.keys()}
for arm in list(arms_dict.keys()):
    for i_night in nights_range:

        # Get the lines for the analysis
        complete_df = logs_dict[arm][i_night]
        idcs_lines = complete_df.index.isin(ext_lines[arm])
        slice_df = complete_df.loc[idcs_lines]

        # Compute the extinction:
        redcorr_dict[arm][i_night] = compute_cHbeta(slice_df, reddening_curve=red_curve, R_v=R_v, ref_wave='H1_4861A_N1',
                                                    compMode='gauss')

# Save extinction values to the configuration file
for arm in list(arms_dict.keys()):
    for i_night in nights_range:
        cHbeta, cHbeta_err = redcorr_dict[arm][i_night]['cHbeta'], redcorr_dict[arm][i_night]['cHbeta_err']
        param_dict = {f'cHbeta_{arm}_{i_night}': [np.round(cHbeta, 3), np.round(cHbeta_err, 3)]}
        safe_cfg(cfg_file, param_dict, section_name='MIKE_extinction_narrow_profile')


# Make the plots
STANDARD_PLOT = {'figure.figsize': (14, 7), 'axes.titlesize': 12, 'axes.labelsize': 14,
                 'legend.fontsize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10}
rcParams.update(STANDARD_PLOT)

axes_dict = {'xlabel': r'$f_{\lambda} - f_{H\beta}$',
             'ylabel': r'$ \left(\frac{I_{\lambda}}{I_{\H\beta}}\right)_{Theo} - \left(\frac{F_{\lambda}}{F_{\H\beta}}\right)_{Obs}$',
             'title': f'Logaritmic extinction coefficient calculation'}

for arm in list(arms_dict.keys()):

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.subplots_adjust(bottom=0.320, top=0.915)
    cmap = cm.get_cmap()

    for i_night in nights_range:

        night_color = cmap(nights_range.index(i_night)/len(nights_range))

        redcorr_night = redcorr_dict[arm][i_night]

        ax.errorbar(redcorr_night['x'], redcorr_night['y'], yerr=redcorr_night['y_err'], fmt='o', color=night_color)

        cHbeta_label = f'{arm.capitalize()} arm, night {i_night}'
        label = r'{}: $c(H\beta)$ = ${:.2f}\pm{:.2f}$'.format(cHbeta_label, redcorr_night['cHbeta'], redcorr_night['cHbeta_err'])

        cHbeta_trndline = redcorr_night['cHbeta'] * redcorr_night['x'] + redcorr_night['intercept']
        ax.plot(redcorr_night['x'], cHbeta_trndline, linestyle='--', label=label, color=night_color)

    ax.update(axes_dict)

    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.7])

    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25),
              fancybox=True, shadow=True)

    # plt.savefig(plot_address, dpi=200)
    plt.show()



