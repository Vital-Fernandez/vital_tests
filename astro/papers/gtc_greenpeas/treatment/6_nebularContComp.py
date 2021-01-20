import numpy as np
import pandas as pd
from pathlib import Path
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.gasContinuum_functions import NebularContinua
import matplotlib.pyplot as plt

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

ext = '_BR'

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

    # Physical parameters
    Te_low = results_dict['Initial_values']['T_low']
    cHbeta = results_dict['Initial_values']['cHbeta_BR_Hbeta_Hgamma_Hdelta']
    HeII_HII, HeIII_HeII = results_dict['Initial_values']['HeII_HII'], results_dict['Initial_values']['HeIII_HII']

    # Load spectrum
    print(f'\n-- Treating: {obj}{ext}.fits')
    wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
    flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
    lm = sr.LineMesurer(wave, flux, redshift=z_array[i], crop_waves=(wmin_array[i], wmax_array[i]))

    # Extinction correction
    rc = pn.RedCorr(R_V=3.4, law='G03 LMC', cHbeta=cHbeta[0])
    red_corr = rc.getCorr(lm.wave)
    int_spec = lm.flux * red_corr

    # Compute nebular continuum
    nebCalc = NebularContinua()
    Hbeta_flux = linesDF.loc['H1_4861A'].intg_flux
    Halpha_redCorr = rc.getCorr(6563)
    Halpha_int = Hbeta_flux * Halpha_redCorr * 2.86
    neb_int = nebCalc.flux_spectrum(lm.wave, Te_low, Halpha_int, HeII_HII, HeIII_HeII)

    # Save object spectrum without nebular component
    flux_noNeb = ((int_spec - neb_int) / red_corr)
    flux_neb = (neb_int/red_corr)
    np.savetxt(nebFluxNoNebCompFile, np.transpose(np.array([lm.wave, flux_noNeb])), fmt="%7.1f %10.4e")
    np.savetxt(nebCompFile, np.transpose(np.array([lm.wave, flux_neb])), fmt="%7.1f %10.4e")

    # Plot spectra components
    labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                  'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                  'title': f'Galaxy {obj}{ext} nebular continuum calculation'}

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(lm.wave, lm.flux, label='Object flux  spectrum')
    ax.plot(lm.wave, int_spec, label='Object intensity spectrum')
    ax.plot(lm.wave, neb_int, label='Nebular intensity spectrum')
    ax.plot(lm.wave, flux_noNeb, label='Object flux no nebular component')
    ax.update(labelsDict)
    ax.legend()
    ax.set_yscale('log')
    # plt.savefig(nebPlotFile, bbox_inches='tight')
    plt.show()
