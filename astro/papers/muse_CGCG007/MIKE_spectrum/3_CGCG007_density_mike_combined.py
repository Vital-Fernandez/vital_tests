import numpy as np
import pyneb as pn
import pandas as pd

from pathlib import Path
from matplotlib import rcParams, pyplot as plt, cm
from astro.papers.muse_CGCG007.muse_CGCG007_methods import read_lines_fits, compute_cHbeta, safe_cfg
from lime.tools import label_decomposition
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

        # Apply the extinction correction
        cHbeta = obsCfg['MIKE_extinction_summed_profiles'][f'cHbeta_{arm}_{i_night}'][0]
        rc = pn.RedCorr(R_V=R_v, law=red_curve, cHbeta=cHbeta)
        wave_array, flux_array, err_array = complete_df.wavelength.values, complete_df.gauss_flux.values, complete_df.gauss_err.values
        corr = rc.getCorr(wave_array)
        intensity, intensityErr = flux_array * corr, err_array * corr
        complete_df['gauss_int'], complete_df['gauss_int_err'] = intensity, intensityErr

# Compute the nSII density
T_low = 15000.0
ne_low = 250
S3 = pn.Atom('S', 3)
S2 = pn.Atom('S', 2)
O2 = pn.Atom('O', 2)
O3 = pn.Atom('O', 3)
mc_steps = 1000

# Loop through the nights
for arm in ['red']:
    for i_night in nights_range:
        nightDF = logs_dict[arm][i_night]
        S2_6716A_dist = np.random.normal(nightDF.loc['S2_6716A'].gauss_int, nightDF.loc['S2_6716A'].gauss_int_err, mc_steps)
        S2_6731A_dist = np.random.normal(nightDF.loc['S2_6731A'].gauss_int, nightDF.loc['S2_6731A'].gauss_int_err, mc_steps)
        nSII = S2.getTemDen(S2_6716A_dist/S2_6731A_dist, tem=T_low, wave1=6716, wave2=6731)
        print(f'Arm {arm}, night {i_night}: ne_SII = {nSII.mean():0.2f}+/-{nSII.std():0.2f} (cm^-3)')

    print()
    for i_night in nights_range:
        nightDF = logs_dict[arm][i_night]
        S3_6312A_dist = np.random.normal(nightDF.loc['S3_6312A'].gauss_int, nightDF.loc['S3_6312A'].gauss_int_err, mc_steps)
        S3_9069A_dist = np.random.normal(nightDF.loc['S3_9069A'].gauss_int, nightDF.loc['S3_9069A'].gauss_int_err, mc_steps)
        TSIII = S3.getTemDen(S3_6312A_dist/S3_9069A_dist, den=ne_low, wave1=6312, wave2=9069)
        print(f'Arm {arm}, night {i_night}: Te_SIII = {TSIII.mean():0.2f}+/-{TSIII.std():0.2f} (K)')

# Loop through the nights
print()
for arm in ['blue']:

    for i_night in nights_range:
        nightDF = logs_dict[arm][i_night]
        O2_3726A_dist = np.random.normal(nightDF.loc['O2_3726A'].gauss_int, nightDF.loc['O2_3726A'].gauss_int_err, mc_steps)
        O2_3729A_dist = np.random.normal(nightDF.loc['O2_3729A'].gauss_int, nightDF.loc['O2_3729A'].gauss_int_err, mc_steps)
        nOII = O2.getTemDen(O2_3726A_dist/O2_3729A_dist, tem=T_low, wave1=3726, wave2=3729)
        print(f'Arm {arm}, night {i_night}: ne_OII = {nOII.mean():0.2f}+/-{nOII.std():0.2f} (cm^-3)')

    print()
    for i_night in nights_range:
        nightDF = logs_dict[arm][i_night]
        O3_4363A_dist = np.random.normal(nightDF.loc['O3_4363A'].gauss_int, nightDF.loc['O3_4363A'].gauss_int_err, mc_steps)
        O3_4959A_dist = np.random.normal(nightDF.loc['O3_4959A'].gauss_int, nightDF.loc['O3_4959A'].gauss_int_err, mc_steps)
        O3_5007A_dist = np.random.normal(nightDF.loc['O3_5007A'].gauss_int, nightDF.loc['O3_5007A'].gauss_int_err, mc_steps)
        TSIII = O3.getTemDen(O3_4363A_dist/O3_5007A_dist, den=ne_low, wave1=4363, wave2=5007)
        print(f'Arm {arm}, night {i_night}: Te_OIII = {TSIII.mean():0.0f}+/-{TSIII.std():0.0f} (K), [OIII] line ratio {np.mean(O3_5007A_dist/O3_4959A_dist):0.2f}')

