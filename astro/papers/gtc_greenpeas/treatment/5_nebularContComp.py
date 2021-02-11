import numpy as np
from pathlib import Path
import pyneb as pn
import src.specsiser as sr
import matplotlib.pyplot as plt
from src.specsiser.physical_model.gasContinuum_functions import NebularContinua
from astro.papers.gtc_greenpeas.common_methods import double_arm_redCorr

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

w_div_array = obsData['sample_data']['w_div']
red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

ext = 'BR'
cycle = 'it1'

# Pyneb objects
H1 = pn.RecAtom('H', 1)

for i, obj in enumerate(objList):

    print(f'Treating: {obj}')

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    fits_file = dataFolder / f'{obj}_{ext}.fits'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'

    # Declare output files
    nebCompFile = objFolder/f'{obj}_{ext}_nebFlux_{cycle}.txt'
    nebPlotFile = objFolder/f'{obj}_{ext}_nebComp_{cycle}.png'

    # Load data
    results_dict = sr.loadConfData(results_file, group_variables=False)

    linesDF = sr.lineslogFile_to_DF(lineLog_file)

    wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
    flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
    lm = sr.LineMesurer(wave, flux, redshift=z_array[i], crop_waves=(wmin_array[i], wmax_array[i]))

    # Physical conditions
    Te_low = results_dict[f'Extinction_{cycle}']['Te_low']
    ne = results_dict[f'Extinction_{cycle}']['ne']
    HeII_HII = results_dict[f'Extinction_{cycle}']['He1r']
    HeIII_HeII = results_dict[f'Extinction_{cycle}']['He2r']

    # Extinction parameters
    cHbeta_label = obsData[obj]['cHbeta_label']
    cHbeta = np.array(results_dict[f'Extinction_{cycle}'][cHbeta_label], dtype=float)
    int_spec, corr_spec = double_arm_redCorr(lm.wave, lm.flux, w_div_array[i], red_law, RV, cHbeta)

    # Calculate Halpha intensity from Hbeta flux
    rc = pn.RedCorr(R_V=RV, law=red_law, cHbeta=cHbeta[0])
    Hbeta_int = linesDF.loc['H1_4861A'].intg_flux * rc.getCorr(4861.0)

    # Compute nebular continuum
    nebCalc = NebularContinua()
    neb_int = nebCalc.flux_spectrum(lm.wave, Te_low, Hbeta_int, HeII_HII, HeIII_HeII)

    # Save nebular flux
    flux_neb = (neb_int/corr_spec)
    np.savetxt(nebCompFile, np.transpose(np.array([lm.wave, flux_neb])), fmt="%7.1f %10.4e")

    # --------------------------- Plot nebular continuum

    labelsDict = {'xlabel': r'Wavelength $(\AA)$',
                  'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
                  'title': f'Galaxy {obj}{ext} nebular continuum calculation {cycle}'}
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(lm.wave, lm.flux, label='Object flux spectrum', color='tab:red')
    ax.plot(lm.wave, flux_neb, label='Nebular flux spectrum', color='tab:red', linestyle='--')
    ax.plot(lm.wave, int_spec, label='Object Intensity spectrum', color='tab:blue')
    ax.plot(lm.wave, neb_int, label='Nebular Intensity spectrum', color='tab:blue', linestyle='--')
    ax.update(labelsDict)
    ax.legend()
    ax.set_yscale('log')
    plt.savefig(nebPlotFile, bbox_inches='tight')
    # plt.show()
