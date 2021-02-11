import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import red_corr_HalphaHbeta_ratio
from uncertainties import ufloat
from uncertainties import ufloat, unumpy
from src.specsiser.physical_model.line_tools import STANDARD_PLOT, STANDARD_AXES
from matplotlib import pyplot as plt, rcParams

objList = ['gp030321', 'gp101157', 'gp121903']
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)

dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']
color_dict = {'_SDSS': 'tab:purple', '_B': 'tab:blue', '_R': 'tab:red'}
label_dict = {'_SDSS': 'SDSS dr16', '_B': 'Blue arm - Osiris', '_R': 'Red arm - Osiris'}

counter = 0
for i, obj in enumerate(objList):

    objFolder = resultsFolder/f'{obj}'
    mainLinesLog = objFolder/f'{obj}_BR_linesLog.txt'

    mainLogDF = sr.lineslogFile_to_DF(mainLinesLog)

    flux_Hbeta, err_Hbeta = mainLogDF.loc['H1_4861A', 'intg_flux'], mainLogDF.loc['H1_4861A', 'intg_err']
    flux_norm = ufloat(flux_Hbeta, err_Hbeta)

    idcs_obsLines = ~mainLogDF.index.str.contains('_b')
    obsLines = mainLogDF.loc[idcs_obsLines].index.values

    linesFlux = mainLogDF.loc[idcs_obsLines, 'intg_flux'].values
    LinesErr = mainLogDF.loc[idcs_obsLines, 'intg_err'].values
    lineFluxArray = unumpy.uarray(linesFlux, LinesErr)
    lineNormArray = lineFluxArray/flux_norm

    ion_array, wave_array, latexLabel_array = sr.label_decomposition(obsLines, combined_dict=obsData['default_line_fitting'])

    # Plot Configuration
    defaultConf = STANDARD_PLOT.copy()
    rcParams.update(defaultConf)
    # plt.locator_params(axis='x', nbins=obsLines.size)

    x_loc = np.arange(obsLines.size)
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xticks(x_loc)
    ax.set_xticklabels(latexLabel_array, rotation=90)
    ax.errorbar(x_loc, unumpy.nominal_values(lineNormArray), yerr=unumpy.std_devs(lineNormArray), fmt='o', color='black',
                label='OSIRIS')

    # Loop through the extensions and put the data in the plot
    for ext in ['_B', '_R']:
        extLinesLog = objFolder / f'{obj}{ext}_linesLog.txt'
        extLogDF = sr.lineslogFile_to_DF(extLinesLog)

        flux_Hbeta_ext, err_Hbeta_ext = extLogDF.loc['H1_4861A', 'intg_flux'], extLogDF.loc['H1_4861A', 'intg_err']
        flux_hbeta_ext = ufloat(flux_Hbeta_ext, err_Hbeta_ext)

        x_arr, y_arr, er_arr = [], [], []
        for i_line, lineLabel in enumerate(obsLines):
            if lineLabel in extLogDF.index:
                flux_line, flux_err = extLogDF.loc[lineLabel, 'intg_flux'], extLogDF.loc[lineLabel, 'intg_err']
                flux_value = ufloat(flux_line, flux_err)
                normFluxValue = flux_value/flux_hbeta_ext

                x_arr.append(i_line)
                y_arr.append(normFluxValue.nominal_value)
                er_arr.append(normFluxValue.std_dev)

        ax.errorbar(x_arr, y_arr, yerr=er_arr, fmt='o', color=color_dict[ext], label=label_dict[ext])

    ax.set_yscale('log')
    ax.update({'ylabel': r'$F_{\lambda}/F_{H\beta}$', 'title': f'Green pea {obj}'})
    ax.legend()
    plt.tight_layout()

    plot_address = objFolder/f'{obj}_line_flux_comparison_log.png'
    # plt.savefig(plot_address, dpi=200, bbox_inches='tight')
    plt.show()

    # Increase counter for obj number
    counter += 1