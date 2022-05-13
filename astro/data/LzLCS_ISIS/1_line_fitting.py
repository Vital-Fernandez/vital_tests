import numpy as np
import lime
from lime.io import load_fits
from pathlib import Path
from shutil import copy as shu_copy


def ordering_mask(log, target_list=['O3_5007A_b', 'O3_4959A_b', 'H1_4861A_b', 'H1_6563A_b']):

    for line in log.index:
        if line not in target_list:
            target_list.append(line)

    log.reindex(target_list, axis=0)

    return log.reindex(target_list, axis=0)


def A_and_K_calculation(log):

    v_5 = log['v_5'].values
    v_10 = log['v_10'].values
    v_90 = log['v_90'].values
    v_95 = log['v_95'].values
    v_med = log['v_med'].values
    FWHM_intg = log['FWHM_intg'].values

    W_80 = v_90 - v_10
    W_90 = v_95 - v_5
    A_factor = ((v_90 - v_med) - (v_med - v_10)) / W_80
    K_factor = W_90 / (1.397 * FWHM_intg)

    peak_waves = log.peak_wave.values
    center_profiles = log.center.values
    trans_waves = log.wavelength.values

    ref_lines = log["profile_label"].str.split('-', expand=True)[0].values
    wave_arr = np.empty(ref_lines.size)
    for i, line_label in enumerate(ref_lines):
        ion_arr, wave_arr[i], latex_array = lime.label_decomposition(line_label, scalar_output=True)

    z_peak = peak_waves/wave_arr - 1
    z_profiles = center_profiles/trans_waves -1
    v_gal = 299792.458 * z_peak
    v_line = 299792.458 * z_profiles
    v_r_fitelp = v_line - v_gal
    v_r_err_fitelp = log.v_r_err.values

    return A_factor, K_factor, W_80, v_r_fitelp, v_r_err_fitelp


conf_file = 'LzLCS_ISIS_cfg.ini'
obsCfg = lime.load_cfg(conf_file, obj_section={'sample_data': 'objName_list'}, def_cfg_sec='default_line_fitting')

dataFolder = Path(obsCfg['data_location']['data_folder'])
treatmentFolder = Path(obsCfg['data_location']['results_folder'])

specNameList = obsCfg['sample_data']['specName_list']
zList = obsCfg['sample_data']['redshift_array']
arm_list = obsCfg['sample_data']['arm_list']
objList = obsCfg['sample_data']['objName_list']
norm_flux = obsCfg['sample_data']['norm_flux']
refMask = '/home/vital/Dropbox/Astrophysics/Data/LzLCS_ISIS/data/reference_mask.txt'

for i, specName in enumerate(specNameList):

    # Loop through the orders
    wave_joined, flux_joined, err_joined = np.array([]), np.array([]), np.array([])
    for arm in arm_list:

        # Input files
        file_spec = dataFolder/objList[i]/f'{specName}_{arm}_f_w_e_flux_nearest.fits'

        # Load the data
        wave, data, hdr = load_fits(file_spec, instrument='ISIS', frame_idx=0)
        flux, err = data[0][0], data[3][0]

        # Crop and join the orders
        wave_joined = np.concatenate([wave_joined, wave])
        flux_joined = np.concatenate([flux_joined, flux])
        err_joined = np.concatenate([err_joined, err])

    # Lime spectrum object
    print(f'- ({i}) {objList[i]}:')
    spec = lime.Spectrum(wave_joined, flux_joined, input_err=err_joined, redshift=zList[i], norm_flux=norm_flux)
    # spec.plot_spectrum(spec_label=objList[i], comp_array=spec.err_flux, frame='rest')

    # Mask file
    objFolder = treatmentFolder/objList[i]
    file_mask = dataFolder/objList[i]/f'{objList[i]}_mask.txt'
    # lime.MaskInspector(file_mask, wave_joined, flux_joined, redshift=zList[i], norm_flux=norm_flux, y_scale='natural',
    #                    n_cols=3, n_rows=2)
    mask = lime.load_lines_log(file_mask)

    # Loop throught the lines
    obj_cfg = obsCfg[f'{objList[i]}_line_fitting']
    for line in mask.index:
        mask_waves = mask.loc[line, 'w1':'w6'].values
        spec.fit_from_wavelengths(line, mask_waves, obj_cfg, fit_method='least_squares')
        spec.display_results(fit_report=False, log_scale=True, output_address=objFolder/f'{line}_gaussian_components.png')
        try:
            spec.plot_line_velocity(output_address=objFolder/f'{line}_velocity_percentiles.png')
        except:
            print(f'This line failed {line}')

    A_array, K_array, w_80_array, v_r_fitelp_arr, v_r_err_fitelp_arr = A_and_K_calculation(spec.log)
    spec.log['A_factor'] = A_array
    spec.log['K_array'] = K_array
    spec.log['w_80'] = w_80_array
    spec.log['v_r_fitelp'] = v_r_fitelp_arr
    spec.log['v_r_err_fitelp'] = v_r_err_fitelp_arr

    # Save line measurements
    lime.save_line_log(spec.log, objFolder/f'{objList[i]}_linesLog.txt')
    lime.save_line_log(spec.log, treatmentFolder/f'ISIS_sample_linesLog.xlsx', ext=objList[i])