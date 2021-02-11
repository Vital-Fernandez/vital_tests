import numpy as np
import pandas as pd
import src.specsiser as sr
import pyneb as pn
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import compute_cHbeta, deredd_fluxes, exitinction_corr_plot, compute_arms_flambda, normalize_flux


conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']
objects_no_chemistry = obsData['file_information']['object_ChemIssues_list']

w_div_array = obsData['sample_data']['w_div']
red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']
H1 = pn.RecAtom('H', 1)

ext = 'SDSS'
cycle = 'it1'

for i, obj in enumerate(objList):

    print(f'- Treating {obj}')

    # Declare input files
    objFolder = resultsFolder/f'{obj}'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'

    # Declare output files
    results_file = objFolder / f'{obj}_BR_measurements.txt'
    plot_address = resultsFolder/f'{obj}'/f'{obj}_{ext}_cHbeta_calculation_{cycle}.png'

    extinction_data = {}
    extinction_data['Te_low'] = 10000.0
    extinction_data['ne'] = 100.0

    # Load the data
    linesDF = sr.lineslogFile_to_DF(lineLog_file)

    # -------------------------------------- Compute cHbeta -------------------------------------
    ext_cor_list, ext_list = [], []

    Te_low, ne = extinction_data['Te_low'], extinction_data['ne']

    extCorrDict = {f'cHbeta_{ext}_Halpha_Hbeta_Hgamma_Hdelta': ['H1_4102A', 'H1_4341A', 'H1_4861A', 'H1_6563A'],
                   f'cHbeta_{ext}_Halpha_Hbeta_Hgamma': ['H1_4341A', 'H1_4861A', 'H1_6563A'],
                   f'cHbeta_{ext}_Hbeta_Hgamma_Hdelta': ['H1_4102A', 'H1_4341A', 'H1_4861A'],
                   f'cHbeta_{ext}_Hbeta_Hgamma': ['H1_4341A', 'H1_4861A']}

    # Fitt reddening coefficient
    reddening_results = compute_cHbeta(extCorrDict, linesDF, red_law, RV, Te_low, ne, compMode='gauss',
                        plot_address=plot_address)

    # Save extinction coefficient for dictionary
    for cHbeta_label, inputLines in extCorrDict.items():
        cHbeta, cHbeta_err = reddening_results[cHbeta_label]['cHbeta'], reddening_results[cHbeta_label]['cHbeta_err']
        extinction_data[cHbeta_label] = (f'{cHbeta:0.3f}', f'{cHbeta_err:0.3f}')
    sr.parseConfDict(results_file, extinction_data, f'Extinction_{ext}_{cycle}', clear_section=True)

    # Plot the data
    # exitinction_corr_plot(obj, ext_cor_list, ext_list)#, plot_save_file=plot_address)

    # -------------------------------------- Apply extinction correction -------------------------------------

    # cHbeta = np.array(extinction_data[cHbeta_label], dtype=float)
    #
    # # Add new columns to dataframe
    # for column in ['obsFlux', 'obsFluxErr', 'f_lambda', 'obsInt', 'obsIntErr']:
    #     linesDF[column] = np.nan
    #
    # # Scale factor for lines in the red arm
    # Halpha_beta_ratio = H1.getEmissivity(tem=Te_low, den=ne, wave=6563) / H1.getEmissivity(tem=Te_low, den=ne, wave=4861)
    #
    # # Distinguish each arm lines
    # idcs_blue = linesDF.wavelength < w_div_array[i]
    # idcs_red = ~idcs_blue
    #
    # # Normalize fluxes blue arm
    # obsFlux, obsFluxErr = normalize_flux(linesDF.loc[idcs_blue], norm_line='H1_4861A')
    # linesDF.loc[idcs_blue, 'obsFlux':'obsFluxErr'] = np.array([obsFlux, obsFluxErr]).T
    #
    # # Normalize fluxes red arm
    # obsFlux, obsFluxErr = normalize_flux(linesDF.loc[idcs_red], norm_line='H1_6563A', scale_factor=Halpha_beta_ratio)
    # linesDF.loc[idcs_red, 'obsFlux':'obsFluxErr'] = np.array([obsFlux, obsFluxErr]).T
    #
    # # Compute reddening law for blue arm
    # f_Hbeta, f_blue = compute_arms_flambda(linesDF.loc[idcs_blue], red_law, RV, ref_line='H1_4861A')
    # linesDF.loc[idcs_blue, 'f_lambda'] = f_blue
    #
    # # Compute reddening law for red arm
    # f_Halpha, f_red = compute_arms_flambda(linesDF.loc[idcs_red], red_law, RV, ref_line='H1_6563A')
    # linesDF.loc[idcs_red, 'f_lambda'] = f_red
    #
    # # Compute line intensities
    # obsInt, obsIntErr = deredd_fluxes(linesDF.obsFlux, linesDF.obsFluxErr, cHbeta[0], cHbeta[1], linesDF.f_lambda)
    # linesDF.loc[:, 'obsInt':'obsIntErr'] = np.array([obsInt, obsIntErr]).T
    #
    # # Save the lines log
    # sr.save_lineslog(linesDF, lineLog_file)
