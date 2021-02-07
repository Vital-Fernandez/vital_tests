import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import compute_arms_flambda, deredd_fluxes, normalize_flux, table_fluxes
import pyneb as pn
from src.specsiser.data_printing import PdfPrinter
import pylatex

conf_file_address = '../../gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']
tables_folder = Path(obsData['file_information']['tables_folder'])

normFlux = obsData['sample_data']['norm_flux']
w_div_array = obsData['sample_data']['w_div']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']
H1 = pn.RecAtom('H', 1)

dict_linesDF = {}
ext = 'BR'
cycle = 'it1'

for i, obj in enumerate(objList):

    print(f'\n- Treating {i} :{obj}_{ext}.fits')

    # Declare input files
    fits_file = dataFolder / f'{obj}_{ext}.fits'
    objFolder = resultsFolder / f'{obj}'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
    results_file = objFolder/f'{obj}_{ext}_measurements.txt'

    # Declare output files
    tables_prefix = objFolder/f'{obj}{ext}_linesTable_{cycle}'

    # Load data
    linesDF = sr.lineslogFile_to_DF(lineLog_file)
    results_dict = sr.loadConfData(results_file, group_variables=False)

    # Physical parameters
    Te = results_dict[f'Extinction_{cycle}']['Te_low']
    ne = results_dict[f'Extinction_{cycle}']['ne']
    cHbeta = results_dict[f'Extinction_{cycle}']['cHbeta_BR_Hbeta_Hgamma_Hdelta']

    # Add new columns to dataframe
    for column in ['obsFlux', 'obsFluxErr', 'f_lambda', 'obsInt', 'obsIntErr']:
        linesDF[column] = np.nan

    # Scale factor for lines in the red arm
    Halpha_beta_ratio = H1.getEmissivity(tem=Te, den=ne, wave=6563) / H1.getEmissivity(tem=Te, den=ne, wave=4861)

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