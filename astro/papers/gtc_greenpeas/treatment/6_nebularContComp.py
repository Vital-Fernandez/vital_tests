import numpy as np
import pandas as pd
from pathlib import Path
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.gasContinuum_functions import NebularContinua
import matplotlib.pyplot as plt

objList = ['gp030321', 'gp101157', 'gp121903']
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList=objList, group_variables=False)

fileList = obsData['file_information']['files_list']
dataFolder = Path(obsData['file_information']['data_folder'])
outputFolder = dataFolder/'flux_analysis'

objList_B = obsData['file_information']['objectB_list']
fileList_B = obsData['file_information']['filesB_list']
objList_R = obsData['file_information']['objectR_list']
fileList_R = obsData['file_information']['filesR_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])


counter = 0
for i, obj in enumerate(objList):

    z = z_objs[i]
    cHbeta = obsData[obj]['cHbeta']
    Te_low = obsData[obj]['Te_low']
    HeII_HII, HeIII_HeII = obsData[obj]['HeII_HII'], obsData[obj]['HeII_HII']

    for ext in ('_BR', '_B'):

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        lineLog_file = outputFolder/f'{obj}{ext}_linesLog.txt'
        nebFluxNoNebCompFile = outputFolder/f'{obj}{ext}_obs_RemoveNebularComp.txt'
        nebCompFile = outputFolder/f'{obj}{ext}_NebFlux.txt'
        nebPlotFile = outputFolder/f'{obj}{ext}_nebComp.png'

        # Set wavelength and flux
        print(f'\n-- Treating {counter} :{obj}{ext}.fits')
        wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        wave_rest = wave / (1 + z)
        idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])
        flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array

        # Load line measurer object
        lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave], linesDF_address=lineLog_file, normFlux=flux_norm)
        nebCalc = NebularContinua()

        # Reddening correction
        rc = pn.RedCorr(R_V=3.4, law='G03 LMC', cHbeta=cHbeta)
        red_corr = rc.getCorr(lm.wave)
        int = lm.flux * red_corr

        # Compute nebular continuum
        Hbeta_flux = lm.linesDF.loc['H1_4861A'].intg_flux
        Halpha_redCorr = rc.getCorr(6563)
        Halpha_int = Hbeta_flux * 2.86 * Halpha_redCorr

        neb_int = nebCalc.flux_spectrum(lm.wave, Te_low, Halpha_int, HeII_HII, HeIII_HeII)

        # Save object spectrum without nebular component
        flux_noNeb = ((int - neb_int) / red_corr) * flux_norm
        flux_neb = (neb_int/red_corr) * flux_norm
        np.savetxt(nebFluxNoNebCompFile, np.transpose(np.array([lm.wave, flux_noNeb])), fmt="%7.1f %10.4e")
        np.savetxt(nebCompFile, np.transpose(np.array([lm.wave, flux_neb])), fmt="%7.1f %10.4e")

        # Plot spectra components
        labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                      'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                      'title': f'Galaxy {obj}{ext} nebular continuum calculation'}

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.plot(lm.wave, lm.flux, label='Object flux  spectrum')
        ax.plot(lm.wave, int, label='Object intensity spectrum')
        ax.plot(lm.wave, neb_int, label='Nebular intensity spectrum')
        ax.plot(lm.wave, flux_noNeb/flux_norm, label='Object flux no nebular component')
        ax.update(labelsDict)
        ax.legend()
        ax.set_yscale('log')
        plt.savefig(nebPlotFile, bbox_inches='tight')

        #ax.plot(lm.wave, int/red_corr, label='Going back', linestyle=':')
        # ax.set_yscale('log')
        # plt.show()