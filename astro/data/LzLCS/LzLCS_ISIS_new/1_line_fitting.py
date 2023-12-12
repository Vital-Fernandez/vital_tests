import numpy as np
import lime
from lime.io import load_fits
from pathlib import Path
from shutil import copy as shu_copy


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


conf_file = 'LzLCS_ISIS_cfg.toml'
obsCfg = lime.load_cfg(conf_file)

dataFolder = Path(obsCfg['data_location']['data_folder'])
treatmentFolder = Path(obsCfg['data_location']['results_folder'])

specNameList = obsCfg['sample_data']['specName_list']
zList = obsCfg['sample_data']['redshift_array']
arm_list = obsCfg['sample_data']['arm_list']
objList = obsCfg['sample_data']['objName_list']
norm_flux = obsCfg['sample_data']['norm_flux']
refMask = '/Users/matiasrodriguez/Desktop/TESIS/VITAL/LzLCS_ISIS_new/data/reference_mask.txt'

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

        # Make the folder to save the data
        objFolder = treatmentFolder/objList[i]
        if not objFolder.exists():
            objFolder.mkdir(parents=True, exist_ok=True)

        # Lime spectrum
        print(f'- ({i}) {objList[i]}:')
        spec = lime.Spectrum(wave_joined, flux_joined, input_err=err_joined, redshift=zList[i], norm_flux=norm_flux)

        # Fit the lines
        bands_df = lime.load_log(dataFolder/objList[i]/f'{objList[i]}_mask.txt')
        spec.fit.frame(bands_df, obsCfg, id_conf_prefix=objList[i], plot_fit=False)
        spec.fit.report()

        # Make the plots
        obj_cfg = obsCfg[f'{objList[i]}_line_fitting']
        for line in bands_df.index:
            spec.plot.bands(line, output_address=objFolder/f'{line}_profile_fitting.png')
            try:
                spec.plot.velocity_profile(line, output_address=objFolder/f'{line}_velocity_percentiles.png')
            except:
                print(f'This line failed {line}')

        A_array, K_array, w_80_array, v_r_fitelp_arr, v_r_err_fitelp_arr = A_and_K_calculation(spec.log)
        spec.log['A_factor'] = A_array
        spec.log['K_array'] = K_array
        spec.log['w_80'] = w_80_array
        spec.log['v_r_fitelp'] = v_r_fitelp_arr
        spec.log['v_r_err_fitelp'] = v_r_err_fitelp_arr

        spec.save_log(objFolder/f'{objList[i]}_linesLog.txt')
        lime.save_log(spec.log, treatmentFolder/f'ISIS_sample_linesLog.xlsx', page=objList[i])


