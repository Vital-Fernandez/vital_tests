import numpy as np
import pandas as pd
import src.specsiser as sr
import pyneb as pn
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import compute_cHbeta
from matplotlib import pyplot as plt, rcParams


def exitinction_corr_plot(objName, corr_dict_list, ext_file_list, ext_color_dict, plot_save_file=None):

    STANDARD_PLOT = {'figure.figsize': (14, 7),
                     'axes.titlesize': 14,
                     'axes.labelsize': 18,
                     'legend.fontsize': 12,
                     'xtick.labelsize': 12,
                     'ytick.labelsize': 12}
    rcParams.update(STANDARD_PLOT)

    fig, ax = plt.subplots(figsize=(8, 4))

    ax.update({'xlabel': r'$f_{\lambda} - f_{H\beta}$',
               'ylabel': r'$ \left(\frac{I_{\lambda}}{I_{\H\beta}}\right)_{Theo} - \left(\frac{F_{\lambda}}{F_{\H\beta}}\right)_{Obs}$',
               'title': f' {objName} logaritmic extinction calculation'})

    for idx_file, corr_lines_dict in enumerate(corr_dict_list):

        ext_file = ext_file_list[idx_file]
        corr_dict, blue_rc = corr_lines_dict['four_lines'], corr_lines_dict['three_lines']

        ax.errorbar(corr_dict['x'], corr_dict['y'], yerr=corr_dict['y_err'], color=color_dict[ext_file], fmt='o')
        all_ylineFit = corr_dict['cHbeta'] * corr_dict['x'] + corr_dict['intercept']
        label = r'$c(H\beta)$ = ${:.2f}\pm{:.2f}$ ' \
                r'($H\alpha$, $H\beta$, $H\gamma$, $H\delta$)'.format(corr_dict['cHbeta'], corr_dict['cHbeta_err'])
        ax.plot(corr_dict['x'], all_ylineFit, color=ext_color_dict[ext_file], label=label, linestyle='--')


        blue_ylineFit = blue_rc['cHbeta'] * blue_rc['x'] + blue_rc['intercept']
        label = r'$c(H\beta)$ = ${:.2f}\pm{:.2f}$ ($H\beta$, $H\gamma$, $H\delta$)'.format(blue_rc['cHbeta'],
                                                                                           blue_rc['cHbeta_err'])
        ax.plot(blue_rc['x'], blue_ylineFit, color=color_dict[ext_file], label=label, linestyle=':')

    ax.legend()
    plt.tight_layout()
    if plot_save_file is not None:
        plt.savefig(plot_address, dpi=200, bbox_inches='tight')
    else:
        plt.show()

    return


conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']
red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']
color_dict = dict(_BR='tab:purple', _B='tab:blue', _SDSS='black')

linesLog_dict = {}
for i, obj in enumerate(objList):

    # Object files
    plot_address = resultsFolder/f'{obj}'/f'{obj}_cHbeta_calculation.png'

    fit_conf = obsData[f'{obj}_line_fitting']
    ext_cor_list, ext_list = [], []

    for ext in ['_BR',  '_B']:

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objFolder = resultsFolder/f'{obj}'
        lineLog_file = objFolder/f'{obj}{ext}_linesLog.txt'
        results_file = objFolder/f'{obj}{ext}_measurements.txt'

        # Create file to safe the measurements
        sr.parseConfDict(results_file, {'SDSS_website': obsData[obj]['SDSS_website']}, 'Metadata')

        # Dictionary to store measurements
        resultsDict = dict(T_low=10000.0, n_e=100.0, HeII_HII=0.01, HeIII_HII=0.001)

        # Check if file exists
        if fits_file.is_file():

            print(f'{obj}{ext}')
            ext_list.append(ext)

            # Load spectrum
            linesDF = sr.lineslogFile_to_DF(lineLog_file)
            Te, ne = resultsDict['T_low'], resultsDict['n_e']

            # Load line measurer object
            if 'H1_4861A' in linesDF.index:

                # Main lines
                all_lines = ['H1_4102A', 'H1_4341A', 'H1_4861A', 'H1_6563A']
                comb_rc = compute_cHbeta(all_lines, linesDF, red_law, RV, Te, ne)

                # Blue arm lines
                blue_arm_lines = ['H1_4102A', 'H1_4341A', 'H1_4861A']
                blue_rc = compute_cHbeta(blue_arm_lines, linesDF, red_law, RV, Te, ne)

                # Store dicts for plotting
                ext_cor_list.append({'four_lines': comb_rc, 'three_lines': blue_rc})

                # Save extinction coefficient for dictionary
                four_lines_ext = (f'{comb_rc["cHbeta"]:0.3f}', f'{comb_rc["cHbeta_err"]:0.3f}')
                resultsDict[f'cHbeta{ext}_Halpha_Hbeta_Hgamma_Hdelta'] = four_lines_ext

                three_lines_ext = (f'{blue_rc["cHbeta"]:0.3f}', f'{blue_rc["cHbeta_err"]:0.3f}')
                resultsDict[f'cHbeta{ext}_Hbeta_Hgamma_Hdelta'] = three_lines_ext

        # Save dictionary with the measurements
        sr.parseConfDict(results_file, resultsDict, 'Initial_values')

    # Plot the data
    exitinction_corr_plot(obj, ext_cor_list, ext_list, color_dict, plot_save_file=None)
