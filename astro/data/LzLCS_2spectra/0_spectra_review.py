import numpy as np
import lime
from pathlib import Path
from shutil import copy as shu_copy
from astropy.io import fits
import os

from matplotlib import pyplot as plt


def open_XSHOOTER_fits(file_address):

    # Open the fits file
    with fits.open(file_address) as hdul:
        data, hdr = hdul[0].data, hdul[0].header

    # for entry, value in hdr.items():
    #     print(entry, value)

    # Reconstruct the wavelength array
    w_min = hdr['CRVAL1']
    dw = hdr['CDELT1']  # dw (Wavelength interval per pixel)
    pixels = hdr['NAXIS1']  # nw number of output pixels
    w_max = w_min + dw * pixels
    wave_arr = np.linspace(w_min, w_max, pixels, endpoint=False) * 10

    return wave_arr, data


def sliceUp(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


conf_file = 'LzLCS_2spectra.ini'
obsCfg = lime.load_cfg(conf_file)

dataFolder = Path(obsCfg['data_location']['data_folder'])
results_fonder = Path(obsCfg['data_location']['results_folder'])
refMask = '/home/vital/Dropbox/Astrophysics/Data/LzLCS_2spectra/reference_mask.txt'

specNameList = obsCfg['sample_data']['specName_list']
zList = obsCfg['sample_data']['redshift_array']
norm_flux = obsCfg['sample_data']['norm_flux']

parent_mask = f'{lime._dir_path}/resources/parent_mask.txt'
parent_df = lime.load_lines_log(parent_mask)

# /home/vital/Dropbox/Astrophysics/Data/LzLCS_2spectra/data/spectra/J123519/1dspectrum_j123519_UVB.fits
# /home/vital/Dropbox/Astrophysics/Data/LzLCS_2spectra/data/J123519/1dspectrum_j123519_UVB.fits
for i, obj in enumerate(specNameList):

    if i == 1:

        order_list = obsCfg['sample_data'][f'order_list']

        # Loop through the orders
        wave_joined, flux_joined, err_joined = np.array([]), np.array([]), np.array([])
        for order_label in order_list:

            specFileName = dataFolder/f'{obj}'/f'1dspectrum_{obj.lower()}_{order_label}.fits'
            errFileName = dataFolder/f'{obj}'/f'1dspectrum_{obj.lower()}_{order_label}_sigma.fits'

            # Load the data
            wave, flux = open_XSHOOTER_fits(specFileName)
            wave, err = open_XSHOOTER_fits(errFileName)

            # Crop and join the orders
            wmin, wmax = obsCfg['sample_data'][f'{obj}_{order_label}_limits_array'] * (1 + zList[i])
            idcs_spec = np.searchsorted(wave, (wmin, wmax))
            wave_joined = np.concatenate([wave_joined, wave[idcs_spec[0]:idcs_spec[1]]])
            flux_joined = np.concatenate([flux_joined, flux[idcs_spec[0]:idcs_spec[1]]])
            err_joined = np.concatenate([err_joined, err[idcs_spec[0]:idcs_spec[1]]])

        # Adjust mask to object
        obj_folder = results_fonder / obj
        obj_mask = obj_folder / f'{obj}_mask.txt'

        spec = lime.Spectrum(wave_joined, flux_joined, input_err=err_joined, redshift=zList[i], norm_flux=norm_flux)
        # spec.plot_spectrum(spec.err_flux, spec_label=f'{obj}', frame='rest')

        # Initial mask
        initial_mask = lime.spectral_mask_generator((spec.wave_rest[0], spec.wave_rest[-1]))
        lines_ref = initial_mask.index.values

        spec.inspect_line_mask(initial_mask, output_log_address=dataFolder/f'{obj}_mask.txt')


