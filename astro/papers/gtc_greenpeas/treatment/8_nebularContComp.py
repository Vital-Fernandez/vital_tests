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
objects_no_chemistry = obsData['file_information']['object_ChemIssues_list']

z_array = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']
arm_wave_boundary = obsData['sample_data']['w_div']

ext = '_BR'
cycle = 'c2'
cycle_ref = 'First_cycle'

# Pyneb objects
H1 = pn.RecAtom('H', 1)

for i, obj in enumerate(objList):

    print(f'Treating: {obj}')

    # Declare files location
    fits_file = dataFolder / f'{obj}{ext}.fits'
    objFolder = resultsFolder / f'{obj}'
    lineLog_file = objFolder / f'{obj}{ext}_linesLog.txt'
    results_file = objFolder / f'{obj}{ext}_measurements.txt'
    nebFluxNoNebCompFile = objFolder/f'{obj}{ext}_obs_RemoveNebularComp_{cycle}.txt'
    nebCompFile = objFolder/f'{obj}{ext}_NebFlux_{cycle}.txt'
    nebPlotFile = objFolder/f'{obj}{ext}_nebComp_{cycle}.png'

    # Load the data
    results_dict = sr.loadConfData(results_file, group_variables=False)
    linesDF = sr.lineslogFile_to_DF(lineLog_file)

    # Physical parameters
    if obj not in objects_no_chemistry:
        Te_low = results_dict[f'{cycle_ref}_Electron_parameters']['Te_low'][0]
        ne = results_dict[f'{cycle_ref}_Electron_parameters']['ne'][0]
        HeII_HII = results_dict[f'{cycle_ref}_Ionic_Abundances']['He1r'][0]
        HeIII_HeII = results_dict[f'{cycle_ref}_Ionic_Abundances']['He2_4686A'][0]
    else:
        Te_low = obsData[obj]['Te_low_array'][0]
        ne = obsData[obj]['ne_array'][0]
        HeII_HII = obsData[obj]['He1_array'][0]
        HeIII_HeII = obsData[obj]['He2_4686A_array'][0]

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
    labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                  'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                  'title': f'Galaxy {obj}{ext} nebular continuum calculation {cycle}'}

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




