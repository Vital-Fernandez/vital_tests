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

ext = 'BR'
cycle = 'it2'

for i, obj in enumerate(objList):

    print(f'- Treating {obj}')
    previousCycle = cycle.replace('2', '1')

    # Declare input files
    objFolder = resultsFolder/f'{obj}'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'

    # Declare output files
    results_file = objFolder/f'{obj}_{ext}_measurements.txt'
    plot_address = resultsFolder/f'{obj}'/f'{obj}_{ext}_cHbeta_calculation_{cycle}.png'

    # Get physical conditions for extinction calculation
    cHbeta_label = f'cHbeta_{ext}_Hbeta_Hgamma_Hdelta'

    extinction_data = {}
    if obj not in objects_no_chemistry:
        extinction_data['Te_low'] = results_file[f'{previousCycle}_electron_parameters']['Te_low'][0]
        extinction_data['ne'] = results_file[f'{previousCycle}_electron_parameters']['n_e'][0]
        extinction_data['He1r'] = results_file[f'{previousCycle}_ionic_Abundances']['He1r'][0]
        extinction_data['He2r'] = results_file[f'{previousCycle}_ionic_Abundances']['He2r'][0]
    else:
        extinction_data['Te_low'] = obsData[obj]['Te_low_array'][0]
        extinction_data['ne'] = obsData[obj]['ne_array'][0]
        extinction_data['He1r'] = obsData[obj]['He1_array'][0]
        extinction_data['He2r'] = obsData[obj]['He2_4686A_array'][0]

    # Load the data
    linesDF = sr.lineslogFile_to_DF(lineLog_file)

    # -------------------------------------- Compute cHbeta -------------------------------------
    ext_cor_list, ext_list = [], []

    Te_low, ne = extinction_data['Te_low'], extinction_data['ne']

    # Main lines
    all_lines = ['H1_4102A', 'H1_4341A', 'H1_4861A', 'H1_6563A']
    comb_rc = compute_cHbeta(all_lines, linesDF, red_law, RV, Te_low, ne)

    # Blue arm lines
    blue_arm_lines = ['H1_4102A', 'H1_4341A', 'H1_4861A']
    blue_rc = compute_cHbeta(blue_arm_lines, linesDF, red_law, RV, Te_low, ne)

    # Store dicts for plotting
    ext_cor_list.append({'four_lines': comb_rc, 'three_lines': blue_rc})
    ext_list.append(ext)

    # Save extinction coefficient for dictionary
    four_lines_ext = (f'{comb_rc["cHbeta"]:0.3f}', f'{comb_rc["cHbeta_err"]:0.3f}')
    extinction_data[f'cHbeta_{ext}_Halpha_Hbeta_Hgamma_Hdelta'] = four_lines_ext

    three_lines_ext = (f'{blue_rc["cHbeta"]:0.3f}', f'{blue_rc["cHbeta_err"]:0.3f}')
    extinction_data[f'cHbeta_{ext}_Hbeta_Hgamma_Hdelta'] = three_lines_ext

    # Save dictionary with the measurements
    sr.parseConfDict(results_file, extinction_data, f'Extinction_{cycle}', clear_section=True)

    # Plot the data
    exitinction_corr_plot(obj, ext_cor_list, ext_list, plot_save_file=plot_address)

    # -------------------------------------- Apply extinction correction -------------------------------------

    cHbeta = np.array(extinction_data[cHbeta_label], dtype=float)

    # Add new columns to dataframe
    for column in ['obsFlux', 'obsFluxErr', 'f_lambda', 'obsInt', 'obsIntErr']:
        linesDF[column] = np.nan

    # Scale factor for lines in the red arm
    Halpha_beta_ratio = H1.getEmissivity(tem=Te_low, den=ne, wave=6563) / H1.getEmissivity(tem=Te_low, den=ne, wave=4861)

    # Distinguish each arm lines
    idcs_blue = linesDF.wavelength < w_div_array[i]
    idcs_red = ~idcs_blue

    # Normalize fluxes blue arm
    obsFlux, obsFluxErr = normalize_flux(linesDF.loc[idcs_blue], norm_line='H1_4861A')
    linesDF.loc[idcs_blue, 'obsFlux':'obsFluxErr'] = np.array([obsFlux, obsFluxErr]).T

    # Normalize fluxes red arm
    obsFlux, obsFluxErr = normalize_flux(linesDF.loc[idcs_red], norm_line='H1_6563A', scale_factor=Halpha_beta_ratio)
    linesDF.loc[idcs_red, 'obsFlux':'obsFluxErr'] = np.array([obsFlux, obsFluxErr]).T

    # Compute reddening law for blue arm
    f_Hbeta, f_blue = compute_arms_flambda(linesDF.loc[idcs_blue], red_law, RV, ref_line='H1_4861A')
    linesDF.loc[idcs_blue, 'f_lambda'] = f_blue

    # Compute reddening law for red arm
    f_Halpha, f_red = compute_arms_flambda(linesDF.loc[idcs_red], red_law, RV, ref_line='H1_6563A')
    linesDF.loc[idcs_red, 'f_lambda'] = f_red

    # Compute line intensities
    obsInt, obsIntErr = deredd_fluxes(linesDF.obsFlux, linesDF.obsFluxErr, cHbeta[0], cHbeta[1], linesDF.f_lambda)
    linesDF.loc[:, 'obsInt':'obsIntErr'] = np.array([obsInt, obsIntErr]).T

    # Save the lines log
    sr.save_lineslog(linesDF, lineLog_file)