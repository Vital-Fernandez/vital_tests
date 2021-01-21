import numpy as np
import pandas as pd
from pathlib import Path
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.gasContinuum_functions import NebularContinua
import matplotlib.pyplot as plt
from astro.papers.gtc_greenpeas.common_methods import double_arm_redCorr, compute_spectrum_flambda, deredd_fluxes, normalize_flux, table_fluxes


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

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']
arm_wave_boundary = obsData['sample_data']['w_div']

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
    HeII_HII, HeIII_HeII = results_dict['Initial_values']['HeII_HII'], results_dict['Initial_values']['HeIII_HII']

    # Load spectrum
    print(f'\n-- Treating: {obj}{ext}.fits')
    wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
    flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
    lm = sr.LineMesurer(wave, flux, redshift=z_array[i], crop_waves=(wmin_array[i], wmax_array[i]))

    # Spectrum extinction correction
    cHbeta = results_dict['Initial_values']['cHbeta_BR_Hbeta_Hgamma_Hdelta']
    int_spec, corr_spec = double_arm_redCorr(lm.wave, lm.flux, arm_wave_boundary[i], red_law, RV, cHbeta)

    # Calculate Halpha intensity from Hbeta flux
    rc = pn.RedCorr(R_V=RV, law=red_law, cHbeta=cHbeta[0])
    Hbeta_int = linesDF.loc['H1_4861A'].intg_flux * rc.getCorr(4861.0)
    emis_AlphaBetaRatio = H1.getEmissivity(tem=Te_low, den=ne, wave=6563)/H1.getEmissivity(tem=Te_low, den=ne, wave=4861)
    Halpha_int = Hbeta_int * emis_AlphaBetaRatio

    # Compute nebular continuum
    nebCalc = NebularContinua()
    neb_int = nebCalc.flux_spectrum(lm.wave, Te_low, Halpha_int, HeII_HII, HeIII_HeII)

    # Save object spectrum without nebular component
    flux_noNeb = ((int_spec - neb_int) / corr_spec)
    flux_neb = (neb_int/corr_spec)
    np.savetxt(nebFluxNoNebCompFile, np.transpose(np.array([lm.wave, flux_noNeb])), fmt="%7.1f %10.4e")
    np.savetxt(nebCompFile, np.transpose(np.array([lm.wave, flux_neb])), fmt="%7.1f %10.4e")

    # Plot spectra components
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(lm.wave, lm.flux, label='Object flux  spectrum')
    ax.plot(lm.wave, int_spec, label='Object intensity spectrum')
    ax.plot(lm.wave, neb_int, label='Nebular intensity spectrum')
    ax.plot(lm.wave, flux_noNeb, label='Object Flux no nebular component', linestyle='--')
    # ax.plot(lm.wave, flux_neb, label='Nebular Flux spectrum backwards', linestyle='--')
    ax.update(labelsDict)
    ax.legend()
    ax.set_yscale('log')
    plt.savefig(nebPlotFile, bbox_inches='tight')
    # plt.show()
















    # # Hbeta red correction
    # Hbeta_flux = linesDF.loc['H1_4861A'].intg_flux
    # f_spec_Hbeta, f_Hbeta = compute_spectrum_flambda(lm.wave, red_law, RV, ref_line='H1_4861A')
    # scale_Hbeta = 1/Hbeta_flux
    # int_spec_Hbeta, int_spec_mine_err = deredd_fluxes(lm.flux * scale_Hbeta, np.zeros(lm.flux.size), cHbeta[0], cHbeta[1], f_spec_Hbeta)
    #
    # # Halpha red correction
    # Halpha_flux = linesDF.loc['H1_6563A'].intg_flux
    # f_spec_Halpha, f_Halpha = compute_spectrum_flambda(lm.wave, red_law, RV, ref_line='H1_6563A')
    # scale_Halpha = 1/Hbeta_flux * emis_AlphaBetaRatio
    # int_spec_Halpha, int_spec_mine_err = deredd_fluxes(lm.flux * scale_Halpha, np.zeros(lm.flux.size), cHbeta[0], cHbeta[1], f_spec_Halpha)
    #
    # fig, ax = plt.subplots(figsize=(12, 8))
    # ax.plot(lm.wave, lm.flux, label='Flux  spectrum')
    # ax.plot(lm.wave, int_spec, label='Intensity spectrum')
    # ax.plot(lm.wave, int_spec_Hbeta * np.power(10, cHbeta[0]) / scale_Hbeta, label='Intensity spectrum HBeta corr', linestyle=':')
    # ax.plot(lm.wave, int_spec_Halpha * (10 ** (0.4 * rc.E_BV * rc.X(6563))) / scale_Halpha, label='Intensity spectrum Halpha corr', linestyle=':')
    # ax.update(labelsDict)
    # ax.legend()
    # ax.set_yscale('log')
    # plt.show()
