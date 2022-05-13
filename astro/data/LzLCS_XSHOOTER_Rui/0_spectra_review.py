import numpy as np
import lime
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
    dw = hdr['CDELT1']  # dw (Wavelength interval per pixel)
    pixels = hdr['NAXIS1']  # nw number of output pixels
    w_max = w_min + dw * pixels
    wave_arr = np.linspace(w_min, w_max, pixels, endpoint=False) * 10

    return wave_arr, data


def sliceUp(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


conf_file = 'LzLCS_XSHOOTER_Rui_cfg.ini'
obsCfg = lime.load_cfg(conf_file)

dataFolder = Path(obsCfg['data_location']['data_folder'])
results_fonder = Path(obsCfg['data_location']['results_folder'])
refMask = '/home/vital/Dropbox/Astrophysics/Data/LzLCS_XSHOOTER_Rui/osiris_mask.txt'

specNameList = obsCfg['sample_data']['specName_list']
zList = obsCfg['sample_data']['redshift_array']
norm_flux = obsCfg['sample_data']['norm_flux']

for i, obj in enumerate(specNameList):

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

    # if not os.path.exists(obj_folder):
    #     os.makedirs(obj_folder)
    # shu_copy(refMask, obj_mask)
    #
    # # # Find lines
    # mask = lime.load_lines_log(obj_mask)
    # peaks_table, matched_masks_DF = spec.match_line_mask(mask, obsCfg['sample_data']['noiseRegion_array'],
    #                                                      detect_threshold=3, emis_threshold=(3, 3),
    #                                                      abs_threshold=(3, 3), width_tol=10)
    #
    # for line in ['S2_6716A_b']:
    #     matched_masks_DF.loc[line, 'w1':'w6'] = mask.loc[line, 'w1':'w6']
    # matched_masks_DF.sort_values(['w3'], ascending=[True], inplace=True)
    #
    # lime.save_line_log(matched_masks_DF, obj_mask)


    # matched_masks_DF.index = matched_masks_DF.index + '_b'
    # idcs_mult = matched_masks_DF.index.isin([H1_3970A, H1_4102A, ])
    # print(matched_masks_DF)

    mask = lime.load_lines_log(obj_mask)
    spec.plot_spectrum(match_log=mask, spec_label=f'{obj} spectrum', log_scale=False, frame='rest')

    for lineIntervals in sliceUp(mask.index.values, 4):
        lime.MaskInspector(obj_mask, wave_joined, flux_joined, redshift=zList[i], norm_flux=norm_flux,
                           n_cols=2, lines_interval=lineIntervals, y_scale='natural')
