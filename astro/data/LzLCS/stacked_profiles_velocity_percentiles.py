import numpy as np

import pandas
import numpy
import lime
from lime import model
from matplotlib import pyplot as plt
from scipy import stats, optimize
from lime.tools import compute_FWHM0
from pathlib import Path


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
    for i, line_label in enumerate(['O3_5007A']):
        ion_arr, wave_arr[i], latex_array = lime.label_decomposition(line_label, scalar_output=True)
    print(wave_arr)

    z_peak = peak_waves/wave_arr - 1
    z_profiles = center_profiles/trans_waves -1
    v_gal = 299792.458 * z_peak
    v_line = 299792.458 * z_profiles
    v_r_fitelp = v_line - v_gal
    v_r_err_fitelp = log.v_r_err.values

    return A_factor, K_factor, W_80, v_r_fitelp, v_r_err_fitelp


c_KMpS = 299792.458  # Speed of light in Km/s (https://en.wikipedia.org/wiki/Speed_of_light)

folder = Path('/home/vital/Dropbox/Astrophysics/Data/LzLCS_stacks/')
output_folder = folder/'velocity_percentiles'

file_list = ['Non-LCE_stack_normalized_data.csv',
             'Non-LCE_stack_normalized_model.csv',
             'Strong LCE_stack_normalized_data.csv',
             'Strong LCE_stack_normalized_model.csv',
             'Weak LCE_stack_normalized_data.csv',
             'Weak LCE_stack_normalized_model.csv']

for i, file_name in enumerate(file_list):

    vel_array, flux_array, err_array = np.loadtxt(f'{folder}/{file_name}', delimiter=',', unpack=True, skiprows=1)

    # Convert to wavelength
    idx_peak = np.argmax(flux_array)
    peak_wave = 5007.0
    wave_array = (5007.0 / c_KMpS) * vel_array + 5007.0

    spec = lime.Spectrum(wave_array, flux_array, input_err=err_array)
    # spec.plot_spectrum(comp_array=err_array)

    # # Fit the line
    line = 'O3_5007A'
    mask = np.array([4993, 4996, 4997, 5018, 5019, 5022])
    spec.fit_from_wavelengths(line, mask)
    spec.plot_line_velocity(output_address=output_folder/f'{file_name[0:-4]}_vel_percentiles.png')

    A_array, K_array, w_80_array, v_r_fitelp_arr, v_r_err_fitelp_arr = A_and_K_calculation(spec.log)
    spec.log['A_factor'] = A_array
    spec.log['K_array'] = K_array
    spec.log['w_80'] = w_80_array
    spec.log['v_r_fitelp'] = v_r_fitelp_arr
    spec.log['v_r_err_fitelp'] = v_r_err_fitelp_arr

    # Save to log
    ext_name = file_name[0:-4].replace('_stack_normalized', '')
    lime.save_line_log(spec.log, output_folder/f'log_stacked.xlsx', ext=ext_name)

    # fig, ax = plt.subplots(figsize=(12, 12))
    # ax.plot(spec.wave, spec.flux, label='Line profile')
    # for wX in mask:
    #     ax.axvline(x=wX, color='black', linestyle='--')
    # ax.legend()
    # ax.update({'xlabel': 'Velocity', 'ylabel': 'Flux', 'title': f'{file_name}'})
    # plt.show()


# spec = lime.Spectrum(vel, flux, input_err=err)
#
# spec.velocity_percentiles_calculations()
    # # Use lime functions to compute the linear continuum
    # spec_vel = lime.Spectrum(wave_array, flux, input_err=err)
    # idcs_line, idcs_cont = spec_vel.define_masks(vel, [-1050, -950, -500, 500, 950, 1050])
    # contVel, contFlux = vel[idcs_cont], flux[idcs_cont]
    # m_cont, n_cont, r_value, p_value, std_err = stats.linregress(contVel, contFlux)
    # cont_flux = m_cont * vel + n_cont
    #
    # # Peak coordinate
    # idx_peak = np.argmax(flux)
    # pixel_vel = np.diff(vel).mean()
    # pixel_width = 5007.0 * pixel_vel/c_KMpS
    #
    # # FWZI band calculation
    # idx_0 = compute_FWHM0(idx_peak, flux, -1, cont_flux, emission_check=True)
    # idx_f = compute_FWHM0(idx_peak, flux, 1, cont_flux, emission_check=True)
    # w_i, w_f = vel[idx_0], vel[idx_f]
    # lineWave, lineFlux, lineCont = vel[idx_0:idx_f], flux[idx_0:idx_f], lineCont[idx_0:idx_f]
    #
    # percentFluxArray = np.cumsum(lineFlux - lineCont) * pixel_width / self.intg_flux * 100
    #
    # # Cumulative sum flux array
    # percentFluxArray = np.cumsum(line_flux - cont_flux) * pixel_width / self.intg_flux * 100
    #
    # fig, ax = plt.subplots(figsize=(12, 12))
    # ax.plot(vel, flux, label='Line profile')
    # ax.plot(vel, cont_flux, linestyle=':', label='Linear continuum')
    # ax.legend()
    # ax.update({'xlabel': 'Velocity', 'ylabel': 'Flux', 'title': f'{file_name}, min-max = {np.min(vel):0.1f}-{np.max(vel):0.1f} km/s'})
    # plt.show()