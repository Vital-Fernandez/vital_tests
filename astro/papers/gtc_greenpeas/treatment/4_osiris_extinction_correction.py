import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import compute_arms_flambda, deredd_fluxes, normalize_flux, table_fluxes
import pyneb as pn

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
w_div_array = obsData['sample_data']['w_div']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])
red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

# Pyneb objects
H1 = pn.RecAtom('H', 1)

counter = 0
for i, obj in enumerate(objList):

    for ext in ['_BR']:

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objFolder = resultsFolder/f'{obj}'
        lineLog_file = objFolder/f'{obj}{ext}_linesLog.txt'
        tables_prefix = objFolder/f'{obj}{ext}_linesTable'
        results_file = objFolder/f'{obj}{ext}_measurements.txt'

        # Check if file exists
        if fits_file.is_file():

            # Load data
            print(f'\n- Treating {counter}: {obj}{ext}.fits')
            linesDF = sr.lineslogFile_to_DF(lineLog_file)
            results_dict = sr.loadConfData(results_file, group_variables=False)

            # Physical parameters
            Te, ne = results_dict['Initial_values']['T_low'], results_dict['Initial_values']['n_e']
            cHbeta = results_dict['Initial_values'][f'cHbeta{ext}_3lines_obs']

            # Add new columns to dataframe
            for column in ['obsFlux', 'obsFluxErr', 'f_lambda', 'obsInt', 'obsIntErr']:
                linesDF[column] = np.nan

            # Scale factor for lines in the red arm
            H_apha_beta_emisRatio = H1.getEmissivity(tem=Te, den=ne, wave=6563) / H1.getEmissivity(tem=Te, den=ne,
                                                                                                   wave=4861)
            # Distinguish each arm lines
            idcs_blue = linesDF.wavelength < w_div_array[i]
            idcs_red = ~idcs_blue

            # Normalize fluxes blue arm
            obsFlux, obsFluxErr = normalize_flux(linesDF.loc[idcs_blue], norm_line='H1_4861A')
            linesDF.loc[idcs_blue, 'obsFlux':'obsFluxErr'] = np.array([obsFlux, obsFluxErr]).T

            # Normalize fluxes red arm
            obsFlux, obsFluxErr = normalize_flux(linesDF.loc[idcs_red], norm_line='H1_6563A', scale_factor=H_apha_beta_emisRatio)
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

            # Save the line intensities
            idcs_obs = ~linesDF.index.str.contains('_b')
            table_fluxes(linesDF.loc[idcs_obs].index,
                         linesDF.loc[idcs_obs].f_lambda,
                         linesDF.loc[idcs_obs].obsFlux,
                         linesDF.loc[idcs_obs].obsFluxErr,
                         linesDF.loc[idcs_obs].obsInt,
                         linesDF.loc[idcs_obs].obsIntErr,
                         cHbeta[0],
                         cHbeta[1],
                         ref_label='H1_4861A',
                         ref_flux=linesDF.loc['H1_4861A'].intg_flux,
                         ref_err=linesDF.loc['H1_4861A'].intg_err,
                         output_address=tables_prefix)

            # Save the lines log
            sr.save_lineslog(linesDF, lineLog_file)

            # Increase counter for obj number
            counter += 1
