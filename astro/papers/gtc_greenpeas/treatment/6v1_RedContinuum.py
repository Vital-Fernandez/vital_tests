import numpy as np
import pandas as pd
from pathlib import Path
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.gasContinuum_functions import NebularContinua
import matplotlib.pyplot as plt
from astro.papers.gtc_greenpeas.common_methods import compute_spectrum_flambda, deredd_fluxes, double_arm_redCorr


conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']
tables_folder = Path(obsData['file_information']['tables_folder'])
idx_band = int(obsData['file_information']['band_flux'])

z_array = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
arm_wave_boundary = obsData['sample_data']['w_div']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

ext = '_BR'

# Pyneb objects
H1 = pn.RecAtom('H', 1)

for i, obj in enumerate(objList):

    # Declare files location
    fits_file = dataFolder / f'{obj}{ext}.fits'
    objFolder = resultsFolder / f'{obj}'
    lineLog_file = objFolder / f'{obj}{ext}_linesLog.txt'
    results_file = objFolder / f'{obj}{ext}_measurements.txt'
    nebFluxNoNebCompFile = objFolder/f'{obj}{ext}_obs_RemoveNebularComp.txt'
    nebCompFile = objFolder/f'{obj}{ext}_NebFlux.txt'
    nebPlotFile = objFolder/f'{obj}{ext}_nebComp.png'
    results_dict = sr.loadConfData(results_file, group_variables=False)
    linesDF = sr.lineslogFile_to_DF(lineLog_file)

    labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                  'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                  'title': f'Galaxy {obj}{ext} nebular continuum calculation'}

    # Physical parameters
    Te_low, ne = results_dict['Initial_values']['T_low'], results_dict['Initial_values']['n_e']
    cHbeta = results_dict['Initial_values']['cHbeta_BR_Hbeta_Hgamma_Hdelta']
    HeII_HII, HeIII_HeII = results_dict['Initial_values']['HeII_HII'], results_dict['Initial_values']['HeIII_HII']

    # Load spectrum
    print(f'\n-- Treating: {obj}{ext}.fits')
    wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
    flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
    lm = sr.LineMesurer(wave, flux, redshift=z_array[i], crop_waves=(wmin_array[i], wmax_array[i]))

    # Extinction correction
    rc = pn.RedCorr(R_V=RV, law=red_law, cHbeta=cHbeta[0])
    red_corr = rc.getCorr(lm.wave)
    red_corr_Hbeta = rc.getCorr(4861.0)
    int_spec = lm.flux * red_corr
    emis_AlphaBetaRatio = H1.getEmissivity(tem=Te_low, den=ne, wave=6563)/H1.getEmissivity(tem=Te_low, den=ne, wave=4861)

    # Hbeta red correction
    f_spec_Hbeta, f_Hbeta = compute_spectrum_flambda(lm.wave, red_law, RV, ref_line='H1_4861A')
    int_spec_Hbeta, int_spec_mine_err = deredd_fluxes(lm.flux, np.zeros(lm.flux.size), cHbeta[0], cHbeta[1], f_spec_Hbeta)

    # Halpha red correction
    f_spec_Halpha, f_Halpha = compute_spectrum_flambda(lm.wave, red_law, RV, ref_line='H1_6563A')
    int_spec_Halpha, int_spec_mine_err = deredd_fluxes(lm.flux, np.zeros(lm.flux.size), cHbeta[0], cHbeta[1], f_spec_Halpha)

    # Mix correction
    idcs_waveBlue = lm.wave < arm_wave_boundary[i]
    idcs_waveRed = lm.wave > arm_wave_boundary[i]
    f_spec_blue, f_Hbeta = compute_spectrum_flambda(lm.wave[idcs_waveBlue], red_law, RV, ref_line='H1_4861A')
    f_spec_red, f_Halpha = compute_spectrum_flambda(lm.wave[idcs_waveRed], red_law, RV, ref_line='H1_6563A')
    int_Blue, int_Blue_err = deredd_fluxes(lm.flux[idcs_waveBlue], np.zeros(lm.flux[idcs_waveBlue].size), cHbeta[0], cHbeta[1], f_spec_blue)
    int_Red, int_Red_err = deredd_fluxes(lm.flux[idcs_waveRed], np.zeros(lm.flux[idcs_waveRed].size), cHbeta[0], cHbeta[1], f_spec_red)

    # Combine the spectra
    int_spec_mix = np.zeros(lm.flux.size)
    int_spec_mix[idcs_waveBlue] = int_Blue * np.power(10, cHbeta[0])
    int_spec_mix[idcs_waveRed] = int_Red * np.power(10, 0.4 * rc.E_BV * rc.X(6563))

    # Apply redening correction according to the arm
    int_array, corr_array = double_arm_redCorr(lm.wave, lm.flux, arm_wave_boundary[i], red_law, RV, cHbeta)

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(lm.wave, lm.flux, label='Flux  spectrum')
    # ax.plot(lm.wave, int_spec, label='Intensity spectrum pyneb')
    ax.plot(lm.wave, int_spec_Hbeta, label='Intensity spectrum relative HBeta', linestyle=':')
    ax.plot(lm.wave, int_spec_Halpha, label='Intensity spectrum relative Halpha', linestyle=':')

    ax.plot(lm.wave, int_spec_Hbeta * np.power(10, cHbeta[0]), label='Intensity spectrum relative HBeta', linestyle=':')
    ax.plot(lm.wave, int_spec_Halpha * np.power(10, 0.4 * rc.E_BV * rc.X(6563)), label='Intensity spectrum relative Halpha', linestyle=':')

    ax.plot(lm.wave, int_spec_mix, label='Joined reddening spectrum', linestyle='--')
    ax.plot(lm.wave, int_array, label='Joined reddening spectrum function', linestyle='--')

    ax.update(labelsDict)
    ax.legend()
    ax.set_yscale('log')
    plt.show()

    # *np.power(10, cHbeta[0])














