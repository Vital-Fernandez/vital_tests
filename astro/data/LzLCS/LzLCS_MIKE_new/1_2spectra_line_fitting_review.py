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


def A_and_K_calculation(log):

    v_5 = log['v_5'].values
    v_10 = log['v_10'].values
    v_50 = log['v_50'].values
    v_90 = log['v_90'].values
    v_95 = log['v_95'].values
    v_med = log['v_med'].values
    FWHM_intg = log['FWHM_intg'].values

    W_80 = v_90 - v_10
    W_90 = v_95 - v_5
    A_factor = ((v_90 - v_50) - (v_50 - v_10)) / W_80
    K_factor = W_90 / (1.397 * FWHM_intg)

    peak_waves = log.peak_wave.values
    center_profiles = log.center.values
    trans_waves = log.wavelength.values

    wave_arr = log.wavelength.to_numpy()

    z_peak = peak_waves/wave_arr - 1
    z_profiles = center_profiles/trans_waves -1
    v_gal = 299792.458 * z_peak
    v_line = 299792.458 * z_profiles
    v_r_fitelp = v_line - v_gal
    v_r_err_fitelp = log.v_r_err.values

    return A_factor, K_factor, W_80, v_r_fitelp, v_r_err_fitelp


def sliceUp(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


conf_file = 'LzLCS_2spectra.toml'
obsCfg = lime.load_cfg(conf_file)#, obj_section={'sample_data': 'specName_list'}, def_cfg_sec='default_line_fitting')

dataFolder = Path(obsCfg['data_location']['data_folder'])
results_fonder = Path(obsCfg['data_location']['results_folder'])

specNameList = obsCfg['sample_data']['specName_list']
zList = obsCfg['sample_data']['redshift_array']
norm_flux = obsCfg['sample_data']['norm_flux']

# target_lines = ['H1_4861.2582A_b', 'O3_4958.8348A_b', 'O3_5006.7664A_b', 'N2_6583.3513A_b', 'S2_6716.3386A_b']
target_lines = ['H1_4861.2582A_b', 'O3_4958.8348A_b', 'O3_5006.7664A_b', 'N2_6583.3513A_b', 'S2_6716.3386A_b']

for i, obj in enumerate(specNameList):

    order_list = obsCfg['sample_data'][f'order_list']
    obj_folder = results_fonder / obj
    mask_file = dataFolder / obj / f'{obj}_mask.txt'

    # Loop through the orders
    wave_joined, flux_joined, err_joined = np.array([]), np.array([]), np.array([])
    for order_label in order_list:

        specFileName = dataFolder/f'{obj}'/f'1dspectrum_{obj.lower()}_{order_label}.fits'
        errFileName = dataFolder/f'{obj}'/f'1dspectrum_{obj.lower()}_{order_label}_sigma.fits'

        # Load the data
        wave, flux = open_XSHOOTER_fits(specFileName)
        wave, err = open_XSHOOTER_fits(errFileName)

        # Crop and join the orders
        wlim_rest = np.array(obsCfg['sample_data'][f'{obj}_{order_label}_limits_array'])
        wmin, wmax = wlim_rest * (1 + zList[i])
        idcs_spec = np.searchsorted(wave, (wmin, wmax))
        wave_joined = np.concatenate([wave_joined, wave[idcs_spec[0]:idcs_spec[1]]])
        flux_joined = np.concatenate([flux_joined, flux[idcs_spec[0]:idcs_spec[1]]])
        err_joined = np.concatenate([err_joined, err[idcs_spec[0]:idcs_spec[1]]])

    # Adjust mask to object
    print(f'\n Fitting {obj}\n')
    spec = lime.Spectrum(wave_joined, flux_joined, input_err=err_joined, redshift=zList[i], norm_flux=norm_flux,)

    spec.fit.frame(mask_file, obsCfg, id_conf_prefix=f'{obj}', progress_output='counter')
    # spec.plot.spectrum(label=f'{obj}', include_fits=True, rest_frame=True)

    # Make the folder to save the data
    objFolder = results_fonder / obj
    if not objFolder.exists():
        objFolder.mkdir(parents=True, exist_ok=True)

    # # Make the plots
    # for line in spec.log.index:
    #     spec.plot.bands(line, output_address=objFolder / f'{line}_profile_fitting.png')
    #     try:
    #         spec.plot.velocity_profile(line, output_address=objFolder / f'{line}_velocity_percentiles.png')
    #     except:
    #         print(f'This line failed {line}')

    A_array, K_array, w_80_array, v_r_fitelp_arr, v_r_err_fitelp_arr = A_and_K_calculation(spec.log)
    spec.log['A_factor'] = A_array
    spec.log['K_array'] = K_array
    spec.log['w_80'] = w_80_array
    spec.log['v_r_fitelp'] = v_r_fitelp_arr
    spec.log['v_r_err_fitelp'] = v_r_err_fitelp_arr

    spec.save_log(objFolder / f'{obj}_linesLog.txt')
    lime.save_log(spec.log, objFolder / f'MIKE_sample_linesLog.xlsx', page=obj)