from pathlib import Path
import numpy as np
import pyneb as pn
import src.specsiser as sr
from src.specsiser.physical_model.gasContinuum_functions import NebularContinua
from src.specsiser.physical_model.chemical_model import TOIII_from_TSIII_relation
import matplotlib.pyplot as plt

# Import the observation data
obsData = sr.loadConfData('../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini', group_variables=False)
# linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
data_folder = Path(obsData['file_information']['data_folder'])
file_list = obsData['file_information']['files_list']
addressList = list(data_folder/file for file in file_list)
flux_norm = obsData['sample_data']['norm_flux']

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # if i == 2:

        # Establish files location
        objName = obsData['file_information']['object_list'][i]
        print(f'-{objName}')

        fitsFolder, fitsFile = file_address.parent, file_address.name
        lineLogFolder = fitsFolder/'flux_analysis'
        lineLogFile, nebPlotFile = fitsFile.replace('.fits', '_linesLog.txt'), fitsFile.replace('.fits', '_nebComp.png')
        nebFluxNoNebCompFile = fitsFile.replace('.fits', '_obs_RemoveNebularComp.txt')
        nebCompFile = fitsFile.replace('.fits', '_NebFlux.txt')

        # Get fits data
        wave, flux, header = sr.import_fits_data(file_address, instrument='OSIRIS')
        z_mean = obsData['sample_data']['z_array'][i]
        wmin_array, wmax_array = obsData['sample_data']['wmin_array'], obsData['sample_data']['wmax_array']

        # Define wave and flux ranges
        wave_rest = wave / (1 + z_mean)
        idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])

        # Load line measurer object
        lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave], lineLogFolder / lineLogFile, normFlux=flux_norm)
        nebCalc = NebularContinua()

        # Apply reddening correction
        rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
        rc.cHbeta = obsData[objName]['cHbeta']
        red_corr = rc.getCorr(lm.wave)
        int = lm.flux * red_corr

        # Compute nebular continuum
        Hbeta_flux = lm.linesDF.loc['H1_4861A'].intg_flux
        Halpha_flux = lm.linesDF.loc['H1_6563A'].gauss_flux
        Halpha_redCorr = rc.getCorr(6563)
        Halpha_int = Hbeta_flux * 2.86 * Halpha_redCorr

        Te_low = obsData[objName]['Te_low']
        HeII_HII, HeIII_HeII = obsData[objName]['HeII_HII'], obsData[objName]['HeII_HII']
        neb_int = nebCalc.flux_spectrum(lm.wave, Te_low, Halpha_int, HeII_HII, HeIII_HeII)

        # Save object spectrum without nebular component
        flux_noNeb = ((int - neb_int) / red_corr) * flux_norm
        flux_neb = (neb_int/red_corr) * flux_norm
        np.savetxt(lineLogFolder / nebFluxNoNebCompFile, np.transpose(np.array([lm.wave, flux_noNeb])), fmt="%7.1f %10.4e")
        np.savetxt(lineLogFolder / nebCompFile, np.transpose(np.array([lm.wave, flux_neb])), fmt="%7.1f %10.4e")

        # Plot spectra components
        labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                      'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                      'title': f'Galaxy {objName} nebular continuum calculation'}

        fig, ax = plt.subplots(figsize=(12, 8))
        ax.plot(lm.wave, lm.flux, label='Object flux  spectrum')
        ax.plot(lm.wave, int, label='Object intensity spectrum')
        ax.plot(lm.wave, neb_int, label='Nebular intensity spectrum')
        ax.plot(lm.wave, flux_noNeb/flux_norm, label='Object flux no nebular component')
        ax.update(labelsDict)
        ax.legend()
        ax.set_yscale('log')
        print('Saving', lineLogFolder/nebPlotFile)
        plt.savefig(lineLogFolder/nebPlotFile, bbox_inches='tight')

        #ax.plot(lm.wave, int/red_corr, label='Going back', linestyle=':')
        # ax.set_yscale('log')
        # plt.show()