import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, red_corr_HalphaHbeta_ratio, default_linelog_types
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
import time
from lmfit.models import PolynomialModel, Model
from src.specsiser.tools.line_fitting import gaussian_model, linear_model, c_KMpS
import fitelp.constants as constants
from fitelp.make_latex_tables import average_velocities_table_to_latex, halpha_regions_table_to_latex
from fitelp.bpt_plotting import bpt_plot
from fitelp.kinematics_calculations import RegionCalculations
from fitelp.fit_line_profiles import plot_profiles
from fitelp.line_profile_info import RegionParameters

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

dict_errs = {}
dict_nan_values = {}

verbose = False

idx_pair = (145, 120)
i = 0
obj = objList[i]

# Data location
cube_address = fitsFolder/fileList[i]
objFolder = resultsFolder/obj
voxelFolder = resultsFolder/obj/'voxel_data'
db_addresss = objFolder / f'{obj}_database.fits'

# Load data
wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')
wave_rest = wave / (1 + z_objs[i])

# Loop throught the line regions
idx_region = 2

# Fitelp dictionary
example_object = RegionParameters(region_name='CGCG07-025',
                                  blue_spec_file=None, #fits file path of the blue spectrum
                                  red_spec_file=None, #fits file path of the red spectrum
                                  blue_spec_error_file=None, #fits file path of the blue spectrum error
                                  red_spec_error_file=None, #fits file path of the red spectrum error
                                  scale_flux=1e20, # Scales the fluxes
                                  center_list={'low': [3649.11211, 3661.84195, 3648.06497]}, # Center values of the Gaussian components for each zone
                                  sigma_list={'low': [1.45, 1.08, 1.14]}, # Sigma values of the Gaussian components for each zone
                                  lin_slope={'low': 0.010116}, # Linear slope values representing the continuum
                                  lin_int={'low': -64.8480}, # Linear intercept values representing the continuum
                                  num_comps={'low': 1}, #Number of Gaussian components for each zone
                                  component_labels=['Narrow 1'], #Labels for each of the gaussian components
                                  component_colors=['b'], #Colour to plot each of the gaussian components
                                  plotting_x_range=[3400, 4000], # xrange of velocities or delta velocities
                                  sigma_instr_blue=4.8,# Instrumental profile on blue arm
                                  sigma_inst_red=4.8,# Instrumental profile on red arm
                                  distance=21.6,  #Distance to object in Mpc
                                  em_lines_for_avg_vel_calc=['OIII-5007A', 'H-beta', 'H-Alpha'], # List of emission-lines to use to calculate the average radial velocity
                                  plot_residuals=True,
                                  show_systemic_velocity=False, # Assumed False if not defined
                                  systemic_velocity=3650)  # Center of most important emission-lines required only if showSystemicVelocity is True


em_line_dict = dict(name=None, plot_color='c', order=1, filter='blue',
                    min_idx=None, max_idx=None,
                    rest_wavelength=None, num_comps=1,
                    amp_list=[15.056061, 4.566674, 5.261243], zone='high', sigma_tsquared=10.39,
                    comp_limits={'a': np.inf, 'c': np.inf, 's': np.inf}, copy_from=None)

name = None
min_idx = None
max_idx = None
rest_wavelength = None

# Voxel mask
region_label = f'region_{idx_region}'
region_mask = fits.getdata(db_addresss, region_label, ver=1)
region_mask = region_mask.astype(bool)

# Lines mask
mask_address = dataFolder / obj / f'{obj}_region{idx_region}_mask.txt'
mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)
user_conf = obsData[f'region{2}_line_fitting']

idx_j, idx_i = idx_pair

flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
flux_err = np.sqrt(cube[:, idx_j, idx_i].var.data) * norm_flux

# Rebin spectrum at ln(Dlambda/lambdaHbeta)
wave_ln = np.log(wave)
InterpolationFlux = interp1d(wave_ln, flux_voxel/norm_flux, kind='slinear', fill_value="extrapolate")
InterpolationErr = interp1d(wave_ln, flux_err/norm_flux, kind='slinear', fill_value="extrapolate")

wave_ln_resample = np.arange(wave_ln.min(), wave_ln.max(), np.diff(wave_ln).min())

# flux_unifSpeed = InterpolationFlux(wave_ln_resample)
# err_unifErr = InterpolationErr(wave_ln_resample)
# wave_unifSpeed = np.power(np.e, wave_ln_resample)
#
#
# # Interpolation = interp1d(wave, flux_voxel, kind='slinear', fill_value="extrapolate")
# # step_lnLmbda = np.diff(wave).mean() / 4861.0
# # wave_lnLambda = 1
# # flux_interpolated = Interpolation(lm.wave)
# # flux_speed = flux_interpolated
#
# fig, ax = plt.subplots()
# ax.step(wave, flux_voxel/norm_flux, label='Observed spectrum in ln(wave)')
# # ax.step(wave_ln_resample, flux_ln_resample, label='Resample in ln(wave)')
# ax.step(wave_unifSpeed, flux_unifSpeed, label='Resample in ln(wave)')
#
# # ax.plot(wave_obs, flux_obs, label='Observed wavelength')
# ax.legend()
# plt.show()

lm = sr.LineMesurer(wave_rest, flux_voxel, input_err=flux_err, redshift=0.0, normFlux=norm_flux)
obsLm = sr.LineMesurer()
zeroLm = sr.LineMesurer()
unifSpeedLm = sr.LineMesurer()

lm.plot_spectrum(specLabel=f'{obj} voxel {idx_j}-{idx_i}', log_scale=False)

# Security check for pixels with nan values:
idcs_nan = np.isnan(lm.flux)
flux_interpolated = None

if idcs_nan.any():
    if region_mask[idx_j, idx_i]:
        Interpolation = interp1d(lm.wave[~idcs_nan], lm.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
        # lm.plot_spectrum_components(continuumFlux=Interpolation(lm.wave))
        flux_interpolated = Interpolation(lm.wave)
        lm.flux = flux_interpolated
        norm_spec = lm.continuum_remover(noise_region)
else:
    norm_spec = lm.continuum_remover(noise_region)

# Identify the emission lines
obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=noise_region, intLineThreshold=3)
maskLinesDF = lm.match_lines(obsLinesTable, mask_df)
idcsObsLines = (maskLinesDF.observation == 'detected')
# lm.plot_spectrum(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF, specLabel=f'Emission line detection')

# Measure the emission lines
lines_list = ['O3_5007A', 'H1_4861A', 'H1_6563A']
for j, lineLabel in enumerate(lines_list):
    if lineLabel in maskLinesDF.index:
        print(f'\n{lineLabel}')

        # # ----------------------- Fitelp fitting
        # wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
        # obsLineWaves = wave_regions * (1 + z_objs[i])
        #
        # example_object.add_em_line(name='OIII-5007A', plot_color='c', order=1, filter='blue', min_idx=2535,
        #                            max_idx=2870, rest_wavelength=5006.84, num_comps=3,
        #                            amp_list=[15.056061, 4.566674, 5.261243], zone='high', sigma_tsquared=10.39,
        #                            comp_limits={'a': np.inf, 'c': np.inf, 's': np.inf}, copy_from=None)
        # example_object.add_em_line(name='H-Alpha', plot_color='y', order=1, filter='blue', min_idx=2800, max_idx=3140,
        #                            rest_wavelength=6562.82, num_comps=3, amp_list=[44.999084, 18.236959, 9.312178],
        #                            zone='low', sigma_tsquared=164.96,
        #                            comp_limits={'a': np.inf, 'c': np.inf, 's': np.inf}, copy_from=None)
        # example_object.add_em_line(name='H-Alpha', plot_color='g', order=4, filter='blue', min_idx=1180,
        #                            max_idx=1510, rest_wavelength=4958.91, num_comps=3,
        #                            amp_list=[5.190979, 1.265695, 0.986356], zone='high', sigma_tsquared=10.39,
        #                            comp_limits={'a': np.inf, 'c': False, 's': False}, copy_from='OIII-5007A')


        # ----------------------- LMFIT fitting
        wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
        lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=user_conf)
        peak_wave = lm.linesDF.loc[lineLabel, "peak_wave"]

        amp = lm.p1[0][0]
        mu = lm.p1[1][0]
        sigma = lm.p1[2][0]

        vr = c_KMpS * (mu - lm.peak_wave) / lm.peak_wave
        sigma_vel = c_KMpS * sigma/lm.peak_wave

        output_message = f'Rest frame: amp ({amp:0.2f}); mu ({mu:0.2f}); sigma ({sigma:0.2f}); vr ({vr:0.2f});  sigma_vel ({sigma_vel:0.2f})'
        print(output_message)

        x_lm, y_lm = lm.fit_output.userkws['x'], lm.fit_output.data
        wave_lm = np.linspace(x_lm[0], x_lm[-1], 500)
        flux_comps = lm.fit_output.eval_components(x=wave_lm)
        flux_lm = flux_comps.get(f'{lineLabel}_cont_', 0.0) + flux_comps[f'{lineLabel}_']

        # ----------------------- Observed frame
        flux_voxel_norm = flux_voxel/norm_flux
        err_voxel_norm = flux_err/norm_flux
        obsLineWaves = wave_regions * (1 + z_objs[i])
        idcsEmis, idcsCont = obsLm.define_masks(wave, flux_voxel_norm, obsLineWaves)

        emisWave, emisFlux = wave[idcsEmis], flux_voxel_norm[idcsEmis]
        contWave, contFlux = wave[idcsCont], flux_voxel_norm[idcsCont]
        obsLm.line_properties(emisWave, emisFlux, contWave, contFlux, bootstrap_size=5000)

        fit_model = Model(linear_model, prefix=f'{lineLabel}_cont_')
        fit_model.set_param_hint(f'{lineLabel}_cont_slope', **{'value': obsLm.m_cont, 'vary':False})
        fit_model.set_param_hint(f'{lineLabel}_cont_intercept', **{'value': obsLm.n_cont, 'vary':False})

        fit_model += Model(gaussian_model, prefix=f'{lineLabel}_')
        fit_model.set_param_hint(f'{lineLabel}_amplitude', value =obsLm.peak_flux - obsLm.cont)
        fit_model.set_param_hint(f'{lineLabel}_center', value = obsLm.peak_wave)
        fit_model.set_param_hint(f'{lineLabel}_sigma', value = 1.0)

        x_array = wave[idcsEmis+idcsCont]
        y_array = flux_voxel_norm[idcsEmis+idcsCont]
        w_array = 1.0 / np.sqrt(err_voxel_norm[idcsEmis+idcsCont])

        fit_params = fit_model.make_params()
        obs_fit_output = fit_model.fit(y_array, fit_params, x=x_array, weights=w_array)

        amp = obs_fit_output.params[f"{lineLabel}_amplitude"].value
        mu = obs_fit_output.params[f"{lineLabel}_center"].value
        sigma = obs_fit_output.params[f"{lineLabel}_sigma"].value

        vr = c_KMpS * (mu - obsLm.peak_wave) / obsLm.peak_wave
        sigma_vel = c_KMpS * sigma/obsLm.peak_wave

        # print(f'intg flux ({obsLm.intg_flux:.3e}); intg err ({obsLm.intg_err:.3e})')
        output_message = f'Observed frame: amp ({amp:0.2f}); mu ({mu:0.2f}); sigma ({sigma:0.2f}); vr ({vr:0.2f});  sigma_vel ({sigma_vel:0.2f})'
        print(output_message)

        x_obs, y_obs = obs_fit_output.userkws['x'], obs_fit_output.data
        wave_obs = np.linspace(x_obs[0], x_obs[-1], 500)
        flux_comps_obs = obs_fit_output.eval_components(x=wave_obs)
        flux_obs = flux_comps_obs.get(f'{lineLabel}_cont_', 0.0) + flux_comps_obs[f'{lineLabel}_']

        # ----------------------- Zero velocity
        flux_voxel_norm = flux_voxel/norm_flux
        err_voxel_norm = flux_err/norm_flux
        obsLineWaves = wave_regions * (1 + z_objs[i])
        idcsEmis, idcsCont = zeroLm.define_masks(wave, flux_voxel_norm, obsLineWaves)

        emisWave, emisFlux = wave[idcsEmis], flux_voxel_norm[idcsEmis]
        contWave, contFlux = wave[idcsCont], flux_voxel_norm[idcsCont]
        idx_peak = np.argmax(emisFlux)
        peak_wave, peak_flux = emisWave[idx_peak], emisFlux[idx_peak]
        emisWave, contWave = c_KMpS*(emisWave - peak_wave)/peak_wave, c_KMpS*(contWave - peak_wave)/peak_wave

        zeroLm.line_properties(emisWave, emisFlux, contWave, contFlux, bootstrap_size=5000)

        zerofit_model = Model(linear_model, prefix=f'{lineLabel}_cont_')
        zerofit_model.set_param_hint(f'{lineLabel}_cont_slope', **{'value': zeroLm.m_cont, 'vary':False})
        zerofit_model.set_param_hint(f'{lineLabel}_cont_intercept', **{'value': zeroLm.n_cont, 'vary':False})
        zerofit_model += Model(gaussian_model, prefix=f'{lineLabel}_')
        zerofit_model.set_param_hint(f'{lineLabel}_amplitude', value=zeroLm.peak_flux - zeroLm.cont)
        zerofit_model.set_param_hint(f'{lineLabel}_center', value=0)
        zerofit_model.set_param_hint(f'{lineLabel}_sigma', value=30.0)

        x_array = c_KMpS * (wave[idcsEmis+idcsCont] - peak_wave)/peak_wave
        y_array = flux_voxel_norm[idcsEmis+idcsCont]
        w_array = 1.0 / np.sqrt(err_voxel_norm[idcsEmis+idcsCont])

        fit_params = zerofit_model.make_params()
        zeroFit_output = zerofit_model.fit(y_array, fit_params, x=x_array, weights=w_array)

        amp = zeroFit_output.params[f"{lineLabel}_amplitude"].value
        mu = zeroFit_output.params[f"{lineLabel}_center"].value
        sigma = zeroFit_output.params[f"{lineLabel}_sigma"].value

        # print(f'intg flux ({obsLm.intg_flux:.3e}); intg err ({obsLm.intg_err:.3e})')
        output_message = f'Velocity frame: amp ({amp:0.2f}); vr ({mu:0.2f}); sigma_vel ({sigma:0.2f})'
        print(output_message)

        x_zero, y_zero = zeroFit_output.userkws['x'], obs_fit_output.data
        wave_zero = np.linspace(x_zero[0], x_zero[-1], 500)
        flux_comps_zero = zeroFit_output.eval_components(x=wave_zero)
        flux_zero = flux_comps_zero.get(f'{lineLabel}_cont_', 0.0) + flux_comps_zero[f'{lineLabel}_']
        x_zero_angs = peak_wave * (1+x_zero/c_KMpS)
        wave_zero_angs = peak_wave * (1+wave_zero/c_KMpS)

        # ----------------------- Uniform velocity frame

        # flux_unifSpeed = InterpolationFlux(wave_ln_resample)
        # err_unifErr = InterpolationErr(wave_ln_resample)
        # wave_unifSpeed = np.power(np.e, wave_ln_resample)

        # obsLineWaves = wave_regions * (1 + z_objs[i])
        # idcsEmis, idcsCont = obsLm.define_masks(wave_unifSpeed, flux_unifSpeed, obsLineWaves)
        #
        # emisWave, emisFlux = wave[idcsEmis], flux_unifSpeed[idcsEmis]
        # contWave, contFlux = wave[idcsCont], flux_unifSpeed[idcsCont]
        # obsLm.line_properties(emisWave, emisFlux, contWave, contFlux, bootstrap_size=5000)
        #
        # unifSpeed_model = Model(linear_model, prefix=f'{lineLabel}_cont_')
        # unifSpeed_model.set_param_hint(f'{lineLabel}_cont_slope', **{'value': obsLm.m_cont, 'vary':False})
        # unifSpeed_model.set_param_hint(f'{lineLabel}_cont_intercept', **{'value': obsLm.n_cont, 'vary':False})
        #
        # unifSpeed_model += Model(gaussian_model, prefix=f'{lineLabel}_')
        # unifSpeed_model.set_param_hint(f'{lineLabel}_amplitude', value = obsLm.peak_flux - obsLm.cont)
        # unifSpeed_model.set_param_hint(f'{lineLabel}_center', value = obsLm.peak_wave)
        # unifSpeed_model.set_param_hint(f'{lineLabel}_sigma', value = 1.0)
        #
        # x_array = wave[idcsEmis+idcsCont]
        # y_array = flux_unifSpeed[idcsEmis+idcsCont]
        # w_array = 1.0 / np.sqrt(err_unifErr[idcsEmis+idcsCont])
        #
        # fit_params = fit_model.make_params()
        # obs_fit_output = fit_model.fit(y_array, fit_params, x=x_array, weights=w_array)
        #
        # amp = obs_fit_output.params[f"{lineLabel}_amplitude"].value
        # mu = obs_fit_output.params[f"{lineLabel}_center"].value
        # sigma = obs_fit_output.params[f"{lineLabel}_sigma"].value
        #
        # vr = c_KMpS * (mu - obsLm.peak_wave)/obsLm.peak_wave
        # sigma_vel = c_KMpS * sigma/obsLm.peak_wave
        #
        # # print(f'intg flux ({obsLm.intg_flux:.3e}); intg err ({obsLm.intg_err:.3e})')
        # output_message = f'Observed frame: amp ({amp:0.2f}); mu ({mu:0.2f}); sigma {sigma:0.2f}; vr ({vr:0.2f});  sigma_vel ({sigma_vel:0.2f})'
        # print(output_message)
        #
        # x_obs, y_obs = obs_fit_output.userkws['x'], obs_fit_output.data
        # wave_obs = np.linspace(x_obs[0], x_obs[-1], 500)
        # flux_comps_obs = obs_fit_output.eval_components(x=wave_obs)
        # flux_obs = flux_comps_obs.get(f'{lineLabel}_cont_', 0.0) + flux_comps_obs[f'{lineLabel}_']
        #
        #
        #
        #
        #
        #
        fig, ax = plt.subplots()
        ax.step(x_obs, y_obs)
        ax.step(x_lm * (1+z_objs[i]), y_lm, linestyle=':')
        ax.plot(wave_obs, flux_obs, label='Observed wavelength')
        ax.plot(wave_lm * (1+z_objs[i]), flux_lm, label='Rest wavelength')

        ax.step(x_zero_angs, y_zero, linestyle='--')
        ax.plot(wave_zero_angs, flux_zero, label='Zero wavelength')
        ax.legend()
        plt.show()


