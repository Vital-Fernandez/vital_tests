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

    ref_lines = log["profile_label"].str.split('-', expand=True)[0].values
    wave_arr = np.empty(ref_lines.size)
    for i, line in enumerate(log.index.values):
        line_label = ref_lines[i]
        if line_label != 'no':
            ion_arr, wave_arr[i], latex_array = lime.label_decomposition(line_label, scalar_output=True)
        else:
            ion_arr, wave_arr[i], latex_array = lime.label_decomposition(line, scalar_output=True)

    z_peak = peak_waves/wave_arr - 1
    z_profiles = center_profiles/trans_waves -1
    v_gal = 299792.458 * z_peak
    v_line = 299792.458 * z_profiles
    v_r_fitelp = v_line - v_gal
    v_r_err_fitelp = log.v_r_err.values

    return A_factor, K_factor, W_80, v_r_fitelp, v_r_err_fitelp


def sliceUp(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))


conf_file = 'LzLCS_2spectra.ini'
obsCfg = lime.load_cfg(conf_file)#, obj_section={'sample_data': 'specName_list'}, def_cfg_sec='default_line_fitting')

dataFolder = Path(obsCfg['data_location']['data_folder'])
results_fonder = Path(obsCfg['data_location']['results_folder'])

specNameList = obsCfg['sample_data']['specName_list']
zList = obsCfg['sample_data']['redshift_array']
norm_flux = obsCfg['sample_data']['norm_flux']

# target_lines = ['H1_4861.2582A_b', 'O3_4958.8348A_b', 'O3_5006.7664A_b', 'N2_6583.3513A_b', 'S2_6716.3386A_b']
target_lines = ['H1_4861.2582A_b', 'O3_4958.8348A_b', 'O3_5006.7664A_b', 'N2_6583.3513A_b', 'S2_6716.3386A_b']

for i, obj in enumerate(specNameList):

    if i == 1:

        order_list = obsCfg['sample_data'][f'order_list']
        obj_folder = results_fonder / obj
        mask_file = dataFolder / f'{obj}_review_mask.txt'

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
        print(f'\n Fitting {obj}\n')
        spec = lime.Spectrum(wave_joined, flux_joined, input_err=err_joined, redshift=zList[i], norm_flux=norm_flux,)
        # spec.plot_spectrum(spec.err_flux, spec_label=f'{obj}', frame='rest')

        # Adjust mask to object
        mask = lime.load_lines_log(mask_file)
        obj_cfg = obsCfg[f'{obj}_review_line_fitting']
        print(obj_cfg)
        for line in target_lines:
            if line in mask.index:
                mask_waves = mask.loc[line, 'w1':'w6'].values
                spec.fit_from_wavelengths(line, mask_waves, obj_cfg, fit_method="least_squares")
                spec.display_results(output_address=obj_folder/f'{line}_gaussian_components.png')
                spec.plot_line_velocity(output_address=obj_folder/f'{line}_kinematics.png')
                spec.display_results(log_scale=True, fit_report=True)

        # A_array, K_array, w_80_array, v_r_fitelp_arr, v_r_err_fitelp_arr = A_and_K_calculation(spec.log)
        # spec.log['A_factor'] = A_array
        # spec.log['K_array'] = K_array
        # spec.log['w_80'] = w_80_array
        # spec.log['v_r_fitelp'] = v_r_fitelp_arr
        # spec.log['v_r_err_fitelp'] = v_r_err_fitelp_arr
        #
        # # Save line measurements
        # lime.save_line_log(spec.log, obj_folder/f'{obj}_linesLog.txt')
        # lime.save_line_log(spec.log, results_fonder/f'2spectra.xlsx', ext=obj)