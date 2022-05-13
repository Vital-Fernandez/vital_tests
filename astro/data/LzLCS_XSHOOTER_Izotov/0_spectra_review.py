import numpy as np
import lime
from lime.io import load_fits
from pathlib import Path
from shutil import copy as shu_copy
from astropy.io import fits
import os

def open_XSHOOTER_fits(file_address):

    # Open the fits file
    with fits.open(file_address) as hdul:
        data, hdr = hdul[0].data, hdul[0].header

    # Reconstruct the wavelength array
    w_min = hdr['CRVAL1']
    dw = hdr['CD1_1']  # dw (Wavelength interval per pixel)
    pixels = hdr['NAXIS1']  # nw number of output pixels
    w_max = w_min + dw * pixels
    wave_arr = np.linspace(w_min, w_max, pixels, endpoint=False)

    return wave_arr, data

conf_file = 'LzLCS_XSHOOTER_Izotov_cfg.ini'
obsCfg = lime.load_cfg(conf_file)

dataFolder = Path(obsCfg['data_location']['data_folder'])
results_fonder = Path(obsCfg['data_location']['results_folder'])

specNameList = obsCfg['sample_data']['specName_list']
zList = obsCfg['sample_data']['redshift_array']
norm_flux = obsCfg['sample_data']['norm_flux']
refMask = '/home/vital/Dropbox/Astrophysics/Data/LzLCS_ISIS/data/reference_mask.txt'

for i, obj in enumerate(specNameList):

    order_list = obsCfg['sample_data'][f'{obj}_order_list']

    # Loop through the orders
    wave_joined, flux_joined, err_joined = np.array([]), np.array([]), np.array([])
    for order_label in order_list:

        specFileName = dataFolder/f'{obj}'/f'f{obj}sum.{order_label}.ms_s.fits'
        errFileName = dataFolder/f'{obj}'/f'f{obj}sum.{order_label}.ms_e.fits'

        # Load the data
        wave, flux = open_XSHOOTER_fits(specFileName)
        wave, err = open_XSHOOTER_fits(errFileName)

        # Crop and join the orders
        wmin, wmax = obsCfg['sample_data'][f'{obj}_{order_label}_array_limits'] * (1 + zList[i])
        idcs_spec = np.searchsorted(wave, (wmin, wmax))
        wave_joined = np.concatenate([wave_joined, wave[idcs_spec[0]:idcs_spec[1]]])
        flux_joined = np.concatenate([flux_joined, flux[idcs_spec[0]:idcs_spec[1]]])

    # LiMe spectrum
    spec = lime.Spectrum(wave_joined, flux_joined, input_err=err, redshift=zList[i], norm_flux=norm_flux)
    # spec.plot_spectrum(spec_label=f'{obj}')

    # Adjust mask to object
    obj_folder = results_fonder/obj
    obj_mask = obj_folder/f'{obj}_mask.txt'

    # if not os.path.exists(obj_folder):
    #     os.makedirs(obj_folder)
    # shu_copy(refMask, obj_mask)

    lime.MaskInspector(obj_mask, wave_joined, flux_joined, redshift=zList[i], norm_flux=norm_flux, n_cols=3)#, y_scale='natural')
