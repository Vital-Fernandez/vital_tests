import numpy as np
import pandas as pd
import src.specsiser as sr
import pyneb as pn
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import red_corr_HalphaHbeta_ratio, compute_cHbeta
from matplotlib import pyplot as plt, rcParams

STANDARD_PLOT = {'figure.figsize': (14, 7),
                 'axes.titlesize': 14,
                 'axes.labelsize': 18,
                 'legend.fontsize': 12,
                 'xtick.labelsize': 12,
                 'ytick.labelsize': 12}
rcParams.update(STANDARD_PLOT)


conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])
red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']
color_dict = dict(_BR='tab:purple', _B='tab:blue', _SDSS='black')

# dictionary with

counter = 0
for i, obj in enumerate(objList):

    z = z_objs[i]
    wmin, wmax = wmin_array[i], wmax_array[i]
    fit_conf = obsData[f'{obj}_line_fitting']

    # Figure for plotting the extinction results
    plot_address = resultsFolder/f'{obj}'/f'{obj}_cHbeta_calculation.png'
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.update({'xlabel': r'$f_{\lambda} - f_{H\beta}$',
               'ylabel': r'$ \left(\frac{I_{\lambda}}{I_{\H\beta}}\right)_{Theo} - \left(\frac{F_{\lambda}}{F_{\H\beta}}\right)_{Obs}$',
               'title': f' {obj} logaritmic extinction calculation'})

    # for ext in ['_BR',  '_B',  '_SDSS']:
    for ext in ['_BR',  '_B']:

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objFolder = resultsFolder/f'{obj}'
        lineLog_file = objFolder/f'{obj}{ext}_linesLog.txt'
        results_file = objFolder/f'{obj}{ext}_measurements.txt'

        # Create file to safe the measurements
        sr.parseConfDict(results_file, {'SDSS_website': obsData[obj]['SDSS_website']}, 'Metadata')

        # Dictionary to store measurements
        resultsDict = dict(T_low=10000.0, n_e=100.0, HeII_HII=100.0, HeIII_HII=0.001)

        if fits_file.is_file():
        # Check if file exists

            print(f'{obj}{ext}')

            # Load spectrum
            linesDF = sr.lineslogFile_to_DF(lineLog_file)
            Te, ne = resultsDict['T_low'], resultsDict['n_e']

            # Load line measurer object
            if 'H1_4861A' in linesDF.index:

                # # All hydrogen lines
                # all_lines = linesDF.loc[linesDF.ion == 'H1'].index.values
                # all_rc = compute_cHbeta(all_lines, linesDF, red_law, RV, Te, ne)
                # ax.errorbar(all_rc['x'], all_rc['y'], yerr=all_rc['y_err'], color='black', fmt='o')
                # ion_array, wavelength_array, latexLabel_array = sr.label_decomposition(all_lines)
                # ax.text(all_rc['x'], all_rc['y'], list(latexLabel_array), fontsize=12, color='black')

                # Main lines
                all_lines = ['H1_4102A', 'H1_4341A', 'H1_4861A', 'H1_6563A']
                comb_rc = compute_cHbeta(all_lines, linesDF, red_law, RV, Te, ne)

                ax.errorbar(comb_rc['x'], comb_rc['y'], yerr=comb_rc['y_err'], color=color_dict[ext], fmt='o')
                all_ylineFit = comb_rc['cHbeta'] * comb_rc['x'] + comb_rc['intercept']
                label = r'$c(H\beta)$ = ${:.2f}\pm{:.2f}$ ' \
                        r'($H\alpha$, $H\beta$, $H\gamma$, $H\delta$)'.format(comb_rc['cHbeta'], comb_rc['cHbeta_err'])
                ax.plot(comb_rc['x'], all_ylineFit, color=color_dict[ext], label=label, linestyle='--')

                # Blue arm lines
                blue_arm_lines = ['H1_4102A', 'H1_4341A', 'H1_4861A']
                blue_rc = compute_cHbeta(blue_arm_lines, linesDF, red_law, RV, Te, ne)

                blue_ylineFit = blue_rc['cHbeta'] * blue_rc['x'] + blue_rc['intercept']
                label = r'$c(H\beta)$ = ${:.2f}\pm{:.2f}$ ($H\beta$, $H\gamma$, $H\delta$)'.format(blue_rc['cHbeta'], blue_rc['cHbeta_err'])
                ax.plot(blue_rc['x'], blue_ylineFit, color=color_dict[ext], label=label, linestyle=':')

                print(f'cHbeta{ext}_4lines_obs = {comb_rc["cHbeta"]:0.3f},{comb_rc["cHbeta_err"]:0.3f}')
                print(f'cHbeta{ext}_3lines_obs = {blue_rc["cHbeta"]:0.3f},{blue_rc["cHbeta_err"]:0.3f}')
                resultsDict[f'cHbeta{ext}_4lines_obs'] = (comb_rc["cHbeta"], comb_rc["cHbeta_err"])
                resultsDict[f'cHbeta{ext}_3lines_obs'] = (blue_rc["cHbeta"], blue_rc["cHbeta_err"])

        # Save dictionary with the measurements
        sr.parseConfDict(results_file, resultsDict, 'Initial_values')


    ax.legend()
    plt.tight_layout()
    plt.savefig(plot_address, dpi=200, bbox_inches='tight')
    # plt.show()
