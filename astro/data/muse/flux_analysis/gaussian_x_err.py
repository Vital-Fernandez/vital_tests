import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from src.specsiser.tools.line_fitting import gaussian_model, linear_model, c_KMpS
from lmfit.models import PolynomialModel, Model
from matplotlib import pyplot as plt

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

idx_j, idx_i = (171, 171)
i = 0
obj = objList[i]
z_galaxy = z_objs[i]
# Data location
cube_address = fitsFolder/fileList[i]
objFolder = resultsFolder/obj
voxelFolder = resultsFolder/obj/'voxel_data'
db_addresss = objFolder / f'{obj}_database.fits'
voxel_spec_file = dataFolder/'voxel_spectrum.txt'

# # Load data
# wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')
# wave_rest = wave / (1 + z_objs[i])
#
# flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
# flux_err = np.sqrt(cube[:, idx_j, idx_i].var.data) * norm_flux
#
# print(voxel_spec_file)
# np.savetxt(voxel_spec_file, np.c_[wave, flux_voxel, flux_err])


# Voxel mask
# Lines mask
mask_address = '/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/CGCG007/CGCG007_region2_mask.txt'
mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)
user_conf = obsData[f'region{3}_line_fitting']

# wave, flux_voxel, flux_err = np.loadtxt(voxel_spec_file, unpack=True)
# lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)
#
# for i, lineLabel in enumerate(['H1_6563A']):
#     wave_regions = mask_df.loc[lineLabel, 'w1':'w6'].values
#     lm.fit_from_wavelengths(lineLabel, wave_regions, user_conf=user_conf)
#     lm.print_results(show_plot=True)

# -------------------------- Another fitting --------------------------------------------

def lmfit_gaussian(x_array, y_array, err_array, x_boundarys):

    # Find indeces for six points in spectrum
    idcsW = np.searchsorted(x_array, x_boundarys)

    # Emission region
    idcsEmis = (x_array[idcsW[2]] <= x_array) & (x_array <= x_array[idcsW[3]])

    # Return left and right continua merged in one array
    idcsCont = (((x_array[idcsW[0]] <= x_array) & (x_array <= x_array[idcsW[1]])) |
                ((x_array[idcsW[4]] <= x_array) & (x_array <= x_array[idcsW[5]]))).squeeze()

    emisWave, emisFlux = x_array[idcsEmis], y_array[idcsEmis]
    contWave, contFlux = x_array[idcsCont], y_array[idcsCont]
    idx_peak = np.argmax(emisFlux)

    fit_model = Model(linear_model, prefix=f'{lineLabel}_cont_')
    fit_model.set_param_hint(f'{lineLabel}_cont_slope', **{'value': 0, 'vary': False})
    fit_model.set_param_hint(f'{lineLabel}_cont_intercept', **{'value': contFlux.mean(), 'vary': False})

    fit_model += Model(gaussian_model, prefix=f'{lineLabel}_')
    fit_model.set_param_hint(f'{lineLabel}_amp', value=emisFlux[idx_peak] - contFlux.mean())
    fit_model.set_param_hint(f'{lineLabel}_center', value=emisWave[idx_peak])
    fit_model.set_param_hint(f'{lineLabel}_sigma', value=1.0)

    x_fit = x_array[idcsEmis + idcsCont]
    y_fit = y_array[idcsEmis + idcsCont]
    w_fit = 1.0 / err_array[idcsEmis + idcsCont]

    fit_params = fit_model.make_params()
    obs_fit_output = fit_model.fit(y_fit, fit_params, x=x_fit, weights=w_fit)

    # amp = obs_fit_output.params[f"{lineLabel}_amp"].value
    mu = obs_fit_output.params[f"{lineLabel}_center"].value
    sigma = obs_fit_output.params[f"{lineLabel}_sigma"].value

    mu_err = obs_fit_output.params[f"{lineLabel}_center"].stderr
    sigma_err = obs_fit_output.params[f"{lineLabel}_sigma"].stderr

    x_obs, y_obs = obs_fit_output.userkws['x'], obs_fit_output.data
    wave_obs = np.linspace(x_obs[0], x_obs[-1], 500)
    flux_comps_obs = obs_fit_output.eval_components(x=wave_obs)
    flux_obs = flux_comps_obs.get(f'{lineLabel}_cont_', 0.0) + flux_comps_obs[f'{lineLabel}_']

    return x_fit, y_fit, wave_obs, flux_obs, mu, sigma, mu_err, sigma_err

lineLabel = 'H1_6563A'
lm = sr.LineMesurer()
wave, flux_voxel, flux_err = np.loadtxt(voxel_spec_file, unpack=True)
flux_voxel, flux_err = flux_voxel/norm_flux, np.sqrt(flux_err/norm_flux)
obsLineWaves = mask_df.loc[lineLabel, 'w1':'w6'].values * (1 + z_objs[i])

MC_size = 500
data_in, curve_out = {}, {}
mu_array, sigma_array = np.zeros(MC_size), np.zeros(MC_size)
mu_err_array, sigma_err_array = np.zeros(MC_size), np.zeros(MC_size)

# # Case all pixels share x_error value
x_err = np.random.normal(loc=0, scale=0.625, size=(MC_size, wave.size))
# x_err = np.random.uniform(low=-0.625, high=0.625, size=MC_size)
for i in range(MC_size):
    x_in, y_in, x_out, y_out, mu, sigma, mu_err, sigma_err = lmfit_gaussian(wave + x_err[i], flux_voxel, flux_err, obsLineWaves)
    data_in[i] = (x_in, y_in)
    curve_out[i] = (x_out, y_out)
    mu_array[i], sigma_array[i] = mu, sigma
    mu_err_array[i], sigma_err_array[i] = mu_err, sigma_err

# Case each pixel has an independent x_error value
x_err = np.random.normal(loc=0, scale=0.625, size=(MC_size, wave.size))
# x_err = np.random.uniform(low=-0.625, high=0.625, size=(MC_size, wave.size))
# for i in range(MC_size):
#     x_in, y_in, x_out, y_out, mu, sigma, mu_err, sigma_err = lmfit_gaussian(wave + x_err[i, :], flux_voxel, flux_err, obsLineWaves)
#     data_in[i] = (x_in, y_in)
#     curve_out[i] = (x_out, y_out)
#     mu_array[i], sigma_array[i] = mu, sigma
#     mu_err_array[i], sigma_err_array[i] = mu_err, sigma_err

# print(f'Guassian center = {mu_array.mean():.2f} +/- {mu_array.std():.2f} (Angstroms) => uncertainty {300000*(6563*(1+z_galaxy)-mu_array).mean()/6563:.2f} km/s')
print(f'Guassian center = {mu_array.mean():.2f} +/- {mu_array.std():.2f} (Angstroms) => uncertainty {300000*mu_array.std()/6563:.2f} km/s')
print(f'Guassian sigma = {sigma_array.mean():.2f} +/- {sigma_array.std():.2f} => uncertainty {300000*sigma_array.std()/6563:.2f} km/s')

fig, ax = plt.subplots()
input_x, input_y = data_in[0]
ax.step(input_x, input_y, where='mid', label=lineLabel)
for i in range(MC_size):
    x_out, f_out = curve_out[i]
    ax.plot(x_out, f_out, alpha=0.1, color='tab:orange')
ax.set_yscale('log')
ax.legend()
plt.show()

# data_in, data_out = {}, {}
# MC_size = 100
# x_err = np.random.uniform(low=-0.625, high=0.625, size=(MC_size, wave.size))
# for i in range(MC_size):
#     x_in, y_in, x_out, y_out = lmfit_gaussian(wave + x_err[i, :], flux_voxel, flux_err, obsLineWaves)
#     data_in[i] = (x_in, y_in)
#     data_out[i] = (x_out, y_out)
#
# fig, ax = plt.subplots()
# input_x, input_y = data_in[0]
# ax.step(input_x, input_y, where='mid', label=r'$H\alpha$')
# for i in range(MC_size):
#     x_out, f_out = data_out[i]
#     ax.plot(x_out, f_out, alpha=0.5)
# ax.set_yscale('log')
# ax.legend()
# plt.show()
