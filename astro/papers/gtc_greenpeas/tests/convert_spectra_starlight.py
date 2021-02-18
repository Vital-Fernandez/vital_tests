import numpy as np
import os
from pathlib import Path
import src.specsiser as sr
from src.specsiser.physical_model.starContinuum_functions import SSPsynthesizer, computeSSP_galaxy_mass
from scipy.interpolate import interp1d
from astro.papers.gtc_greenpeas.common_methods import double_arm_redCorr
import pyneb as pn

if os.name != 'nt':
    conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
    dataFolder = Path('/home/vital/Dropbox/Astrophysics/Papers/gtc_greenpeas/data')
    resultsFolder = Path('/home/vital/Dropbox/Astrophysics/Papers/gtc_greenpeas/treatment')
    starlight_folder = Path('/home/vital/Dropbox/Astrophysics/Tools/Starlight')
else:
    conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
    dataFolder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data')
    resultsFolder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/treatment')
    starlight_folder = Path('D:/Dropbox/Astrophysics/Tools/Starlight')

obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)
objList = obsData['file_information']['object_list']

fileList = obsData['file_information']['files_list']
idx_band = int(obsData['file_information']['band_flux'])

z_array = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
arm_wave_boundary = obsData['sample_data']['w_div']

w_div_array = obsData['sample_data']['w_div']
red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

ext = 'BR'

for i, obj in enumerate(objList):

    print(f'\n-- Treating: {obj}{ext}.fits')

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    fits_file = dataFolder / f'{obj}_{ext}.fits'
    objMask = dataFolder / f'{obj}_{ext}_mask.txt'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_it3.txt'

    # output spectrum
    spectrum_ssp = objFolder / f'{obj}_starlight_spectrum.txt'

    # Load the data
    wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
    flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
    lm = sr.LineMesurer(wave, flux, redshift=z_array[i], crop_waves=(wmin_array[i], wmax_array[i]))
    linesDF = sr.lineslogFile_to_DF(lineLog_file)

    # Add new masks
    linesDF = sr.lineslogFile_to_DF(lineLog_file)
    ini_mask, end_points = obsData[obj]['ini_mask_array'], obsData[obj]['end_mask_array']
    labels = ['cont_mask_' + str(int(x)) for x in ini_mask]
    for j, label_mask in enumerate(labels):
        linesDF.loc[labels[j], ['w3', 'w4']] = ini_mask[j], end_points[j]

    # Compute the mask array
    mask_array = np.zeros(lm.wave.size)

    for mask_label in linesDF.index:
        w_blue, w_red = linesDF.loc[mask_label, ['w3', 'w4']].values
        idcs_mask = (lm.wave > w_blue) & (lm.wave < w_red)
        mask_array[idcs_mask] = 99

    # Compute std error
    idcs_mask = (lm.wave > 4500.0) & (lm.wave < 4550.0)
    err_value = np.std(lm.flux[idcs_mask])
    # error_vector = np.ones(lm.wave.size) * err_value
    error_vector = lm.flux * 0.10

    # Save the spectrum
    array_spectrum = np.transpose(np.array([lm.wave, lm.flux, error_vector, mask_array]))
    np.savetxt(spectrum_ssp, array_spectrum, fmt="%7.1f %10.4e %10.4e %7.0f")

    # Model constants
    c_kms = 300000  # km/s ufloat(74.3, 6.0)
    Huble_Constant = 67.8  # (km/s / Mpc)

    # Compute the distance
    dist_mpc = z_array[i] * (c_kms / Huble_Constant)
    print(f' - {dist_mpc:0.1f}')

