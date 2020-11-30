import numpy as np
import pandas as pd
from pathlib import Path
from src.specsiser.physical_model.starContinuum_functions import SSPsynthesizer
from src.specsiser.physical_model.gasContinuum_functions import NebularContinua
from src.specsiser.inference_model import displaySimulationData
from matplotlib import pyplot as plt, rcParams
import scipy as spy
import theano.tensor as tt
import exoplanet as xo
import pymc3 as pm
import arviz as az

def plot_this(db, ids_db, wave_array, flux_matrix, interp_flux=None):

    fig, ax = plt.subplots(figsize=(9, 5))

    for i_db in db.loc[ids_db].index:
        label = f'{db.loc[i_db].file_name[0:5]}: ' \
                f'log(age) = {np.log10(db.loc[i_db].age_yr)},' \
                f' z_star = {db.loc[i_db].z_star}'
        ax.step(wave_array, flux_matrix[i_db, :], label=label)

    if interp_flux is not None:
        ax.step(wave_array, interp_flux, label='interpolation', linestyle=':')

    ax.set_xlabel(r'$wave (\AA)$')
    ax.set_ylabel(r'$flux norm$')
    ax.legend()
    plt.show()

    return

def generate_region_templates(db, idcs_db, wave_array, flux_array):

    age_range = np.unique(db.loc[idcs_db, 'age_yr'].values)
    z_range = np.unique(db.loc[idcs_db, 'z_star'].values)

    flux_matrix = np.empty((z_range.size, age_range.size, wave_array.size))
    flux_matrix[:, :, :] = np.nan
    for i_z, z in enumerate(z_range):
        for j_age, age in enumerate(age_range):
            i_flux = (db.z_star == z) & (db.age_yr == age)
            flux_matrix[i_z, j_age, :] = flux_array[i_flux]

    return z_range, age_range, flux_matrix



storing_folder = Path('D:\Dropbox\Astrophysics\Papers\gtc_greenpeas\data')

size_dict = {'axes.titlesize': 14, 'axes.labelsize': 14, 'legend.fontsize': 14,
             'xtick.labelsize': 12, 'ytick.labelsize': 12, }
rcParams.update(size_dict)

starCalc = SSPsynthesizer()
nebCalc = NebularContinua()

bases_folder = Path("D:\Dropbox\Astrophysics\Papers\gtc_greenpeas\data\starlight\Bases")
bases_file = Path("D:\Dropbox\Astrophysics\Papers\gtc_greenpeas\data\starlight\Dani_Bases_short.txt")

bases_df_file = storing_folder/'bases_db'
# bases_df, wave_bases, flux_bases = starCalc.import_STARLIGHT_bases(bases_file, bases_folder, crop_waves=None, resam_inter=1, norm_waves=(5100, 5150))
#
# # Storing to files
# bases_df.to_pickle(f'{bases_df_file}_pickle')
# np.save(storing_folder/'waves_bases', wave_bases)
# np.save(storing_folder/'flux_bases', flux_bases)
# bases_df.to_csv(bases_df_file)


# Generate nebular data grid
temp_range = np.linspace(5000, 30000, int((30000-5000)/200)+1)
HeI_range = np.linspace(0.05, 0.20, 51)
HeII_HI = 0.001
# neb_gamma = np.empty((temp_range.size, HeI_range.size, wave_bases.size))
# neb_gamma[:, :, :] = np.nan
# for j, temp in enumerate(temp_range):
#     for i, HeI_HI in enumerate(HeI_range):
#         neb_gamma[j, i, :] = nebCalc.gamma_spectrum(wave_bases, temp, HeI_HI, HeII_HI)
# np.save(storing_folder/'neb_gamma', neb_gamma)


bases_df = pd.read_csv(bases_df_file, index_col=0)
wave_bases = np.load(storing_folder/'waves_bases.npy')
flux_bases = np.load(storing_folder/'flux_bases.npy')
nebular_bases = np.load(storing_folder/'neb_gamma.npy')



# gamma_nebularfb_HI = nebCalc.freebound_gamma(wave, Te, self.HI_fb_dict)
# neb_int = nebCalc.flux_spectrum(lm.wave, Te_low, Halpha_int, HeII_HII, HeIII_HeII)


# starlight_output_file = Path('D:\Dropbox\Astrophysics\Papers\gtc_greenpeas\data\starlight\Output\gp121903_BR.slOutput')
# Input_Wavelength, Input_Flux, Output_Flux, fit_output = starCalc.load_starlight_output(starlight_output_file)
#
# Region 1 templates:
idcs_1 = (bases_df.age_yr < 1.5e7) & (bases_df.z_star > 0.005)
z_range_1, age_range_1, flux_matrix_1 = generate_region_templates(bases_df, idcs_1, wave_bases, flux_bases)
logAge_range_1 = np.log10(age_range_1)
spy_RGridInterp = spy.interpolate.RegularGridInterpolator((z_range_1, logAge_range_1), flux_matrix_1)
exop_interpAxis = xo.interp.RegularGridInterpolator([z_range_1, logAge_range_1], flux_matrix_1)
neb_interpAxis = xo.interp.RegularGridInterpolator([temp_range, HeI_range], nebular_bases)

# Compute nebular component
Temp_true, HeI_II_true, Halpha_True = 12350.0, 0.1055, 500.0
coordB = np.stack(([Temp_true], [HeI_II_true]), axis=-1)
gamma_neb_inter = neb_interpAxis.evaluate(coordB).eval()[0]
flux_neb_inter = nebCalc.zanstra_calibration(wave_bases, Temp_true, Halpha_True, gamma_neb_inter)

# Stellar component
age_true = np.array([6.95])
z_true = np.array([0.0265])
w_true = np.array([1.25])
age_neighbours = [8913000]
z_neighbours = [0.0315321]

synth_spec = np.zeros(wave_bases.size)
spec_components = np.zeros((len(age_true), wave_bases.size))
for i in range(len(age_true)):
    coordB = np.stack(([z_true[i]], [age_true[i]]), axis=-1)
    flux_true_i = exop_interpAxis.evaluate(coordB).eval()
    synth_spec += flux_true_i[0] * w_true[i]
    spec_components[i, :] = flux_true_i[0]

obj_cont = flux_neb_inter + synth_spec
wave_obj, flux_obj_norm, normFlux_coeff_i = starCalc.treat_input_spectrum(wave_bases, obj_cont, norm_waves=(5100, 5150))

# Adding gas and stellar componts
fig, ax = plt.subplots(figsize=(9, 5))
ax.step(wave_obj, synth_spec, label='stellar continuum')
ax.step(wave_obj, flux_neb_inter, label='gas continuum', linestyle=':')
ax.step(wave_obj, flux_neb_inter + synth_spec, label='total continuum', linestyle=':')
ax.step(wave_obj, nebCalc.zanstra_calibration_tt(wave_bases, Temp_true, 500.0, gamma_neb_inter), label='neb zanstra', linestyle=':')
ax.set_xlabel(r'$wave (\AA)$')
ax.set_ylabel(r'$flux norm$')
ax.legend()
plt.show()

# ----- Model with stellar and nebular continuum

# Pymc3 fitter
n_bases = age_true.size
regions_range = np.arange(n_bases)
spec_tensor = tt.zeros(wave_obj.size)
pixelNoise = 0.05 * np.ones(wave_bases.size)

print('Z limits', np.min(z_range_1), np.max(z_range_1))
print('Age limits', np.min(logAge_range_1), np.max(logAge_range_1))

nebular_comp = True
with pm.Model() as model:

    # Define priors
    z_prior = pm.Uniform('z_star', lower=np.min(z_range_1), upper=np.max(z_range_1), shape=n_bases)
    age_prior = pm.Uniform('age', lower=np.min(logAge_range_1), upper=np.max(logAge_range_1), shape=n_bases)
    w_prior = w_true
    Te_prior = pm.Normal('temp', mu=Temp_true, sigma=50) if nebular_comp else Temp_true
    y_plus = pm.Normal('HeI_HI', mu=HeI_II_true, sigma=0.005) if nebular_comp else HeI_II_true

    # Reset continuum container
    if nebular_comp:
        coord_temp = tt.stack([[Te_prior], [y_plus]], axis=-1)
        neb_gamma_t = neb_interpAxis.evaluate(coord_temp)
        spec_tensor = nebCalc.zanstra_calibration_tt(wave_obj, Te_prior, Halpha_True, neb_gamma_t[0])
    else:
        spec_tensor = spec_tensor * 0.0
    # pixelNoise = pm.HalfCauchy('pixelNoise', 0.05)

    # Loop through the stellar components
    for i in regions_range:
        coord_i = tt.stack([[z_prior[i]], [age_prior[i]]], axis=-1)
        spec_i = exop_interpAxis.evaluate(coord_i)
        spec_tensor += w_true[i] * spec_i[0]

    # Likelihood
    pm.Normal('continuum', spec_tensor, pixelNoise, observed=obj_cont)

    # Check simulation statistics
    displaySimulationData(model)

    # Run sampler
    trace = pm.sample(draws=5000, tune=3000, chains=2, cores=1)

print('True values', z_true, age_true, Temp_true, HeI_II_true)

print(az.summary(trace))
az.plot_trace(trace)
plt.show()
az.plot_posterior(trace)
plt.show()
az.plot_forest(trace)
plt.show()





# Compare nebular interpolation
# Temp_true, HeI_II_true = 12350.0, 0.1055
# coordB = np.stack(([Temp_true], [HeI_II_true]), axis=-1)
# flux_neb_inter = neb_interpAxis.evaluate(coordB).eval()[0]
# flux_neb_calc = nebCalc.gamma_spectrum(wave_bases, Temp_true, HeI_II_true, HeII_HI)
# fig, ax = plt.subplots(figsize=(9, 5))
# ax.step(wave_bases, flux_neb_inter, label='flux_neb_inter')
# ax.step(wave_bases, flux_neb_calc, label='flux_neb_calc', linestyle=':')
# ax.set_xlabel(r'$wave (\AA)$')
# ax.set_ylabel(r'$flux norm$')
# ax.legend()
# plt.show()

# # Create synthetic observation
# age_true = np.array([6.95, 6.38])
# z_true = np.array([0.0265, 0.015])
# w_true = np.array([2.5, 3.0])
# age_neighbours = [8913000, 3000000]
# z_neighbours = [0.0315321, 0.019]

# Working
# age_true = np.array([6.95])
# z_true = np.array([0.0265])
# w_true = np.array([1.0])
# age_neighbours = [8913000]
# z_neighbours = [0.0315321]

# age_true = np.array([6.37])
# z_true = np.array([0.009])
# w_true = np.array([1.0])
# age_neighbours = [8913000]
# z_neighbours = [0.0315321]

# synth_spec = np.zeros(wave_bases.size)
# spec_components = np.zeros((len(age_true), wave_bases.size))
# for i in range(len(age_true)):
#     coordB = np.stack(([z_true[i]], [age_true[i]]), axis=-1)
#     flux_true_i = exop_interpAxis.evaluate(coordB).eval()
#     synth_spec += flux_true_i[0] * w_true[i]
#     spec_components[i, :] = flux_true_i[0]
# wave_synth, synth_spec_norm, normFlux_coeff_i = starCalc.treat_input_spectrum(wave_bases, synth_spec, norm_waves=(5100, 5150))


# Compare weights
# fig, ax = plt.subplots(figsize=(9, 5))
# for i in range(len(age_true)):
#     ax.step(wave_synth, spec_components[i, :], label=f' Component {i}')
# ax.step(wave_synth, synth_spec, label='Combined synth')
# ax.step(wave_synth, synth_spec_norm, label='Combined synth norm')
# recons = (spec_components[0] * w_true[0] + spec_components[1] * w_true[1])
# recons_norm = (spec_components[0] * w_true[0]/normFlux_coeff_i + spec_components[1] * w_true[1]/normFlux_coeff_i)
# print(w_true[0]/normFlux_coeff_i + w_true[1]/normFlux_coeff_i, w_true[0]/normFlux_coeff_i, w_true[1]/normFlux_coeff_i)
# ax.step(wave_synth, recons, label='Recons synth')
# ax.step(wave_synth, recons_norm, label='Recons synth norm')
#
# ax.set_xlabel(r'$wave (\AA)$')
# ax.set_ylabel(r'$flux norm$')
# ax.legend()
# plt.show()



# # Pymc3 fitter
# n_bases = age_true.size
# regions_range = np.arange(n_bases)
# spec_tensor = tt.zeros(wave_synth.size)
# pixelNoise = 0.001 * np.ones(wave_bases.size)
#
# print('Z limits', np.min(z_range_1), np.max(z_range_1))
# print('Age limits', np.min(logAge_range_1), np.max(logAge_range_1))
#
# with pm.Model() as model:
#
#     # Reset continuum container
#     spec_tensor = spec_tensor * 0.0
#
#     # Define priors
#     AgeDistBounded = pm.Bound(pm.Normal, lower=np.min(logAge_range_1), upper=np.max(logAge_range_1))
#     MetDistBounded = pm.Bound(pm.Normal, lower=np.min(z_range_1), upper=np.max(z_range_1))
#     z_prior = MetDistBounded('z_star', mu=0.02, sigma=0.01, shape=n_bases)
#     age_prior = AgeDistBounded('age', mu=6.5, sigma=0.05, shape=n_bases)
#     # z_prior = pm.Uniform('z_star', lower=np.min(z_range_1), upper=np.max(z_range_1), shape=n_bases)
#     # age_prior = pm.Uniform('age', lower=np.min(logAge_range_1), upper=np.max(logAge_range_1), shape=n_bases)
#     w_prior = w_true
#     # pixelNoise = pm.HalfCauchy('pixelNoise', 0.05)
#
#     # Loop through the regions
#     for i in regions_range:
#         coord_i = tt.stack([[z_prior[i]], [age_prior[i]]], axis=-1)
#         spec_i = exop_interpAxis.evaluate(coord_i)
#         spec_tensor += w_true[i] * spec_i[0]
#
#     # Likelihood
#     pm.Normal('continuum', spec_tensor, pixelNoise, observed=synth_spec_norm)
#
#     # Check simulation statistics
#     displaySimulationData(model)
#
#     # Run sampler
#     trace = pm.sample(draws=5000, tune=3000, chains=2, cores=1)
#
# print('True values', z_true, age_true)
#
# print(az.summary(trace))
# az.plot_trace(trace)
# plt.show()
# az.plot_posterior(trace)
# plt.show()
# az.plot_forest(trace)
# plt.show()

# # Plot synthetic spectrum with near grid points
# idcs_lib = (bases_df.z_star.isin(z_neighbours)) & (bases_df.age_yr.isin(age_neighbours))
# plot_this(bases_df, idcs_lib, wave_bases, flux_bases, interp_flux=synth_spec_norm)

# # Plot interpolation points
# for z_value in z_range_1:
#     for age_value in age_range_1:
#         print(z_value, age_value)
#         idcs_lib = (bases_df.z_star == z_value) & (bases_df.age_yr == age_value)
#         coordB = np.stack(([z_value], [np.log10(age_value)]), axis=-1)
#         xo_flux_i = exop_interpAxis.evaluate(coordB).eval()
#         plot_this(bases_df, idcs_lib, wave_bases, flux_bases, interp_flux=xo_flux_i[0])

# # Compare scipy and xo interpolation
# z_value, age_value = 0.0075640, 3e6
# coordB = np.stack(([z_value], [age_value]), axis=-1)
# scipy_flux = spy_RGridInterp([[z_value, age_value]])
# xo_flux = exop_interpAxis.evaluate(coordB).eval()
# print('2 Scipy RegularGridInterpolator', scipy_flux)
# print('3 Exoplanet Axis Interpolation', xo_flux)
# idcs_lib = (bases_df.z_star == z_value) & (bases_df.age_yr == age_value)
# plot_this(bases_df, idcs_lib, wave_bases, flux_bases, interp_flux=scipy_flux[0])
# plot_this(bases_df, idcs_lib, wave_bases, flux_bases, interp_flux=xo_flux[0])

0.0190000, 6.749968083509403
0.005, 6.477121254719663
# # Plot by library
# for z_value in [0.0037047, 0.0075640, 0.0190000, 0.0315321]:
#     idcs_lib = (bases_df.z_star == z_value) & (bases_df.age_yr < 1.5e7)
#     plot_this(bases_df, idcs_lib, wave_bases, flux_bases)

# # Plot by library
# for age_value in [1e6, 3e6, 3.9811e+06, 5.6230e+06, 8.9130e+06]:
#     idcs_lib = (bases_df.age_yr == age_value)
#     plot_this(bases_df, idcs_lib, wave_bases, flux_bases)
#
# # Plot by library
# idcs_lib = bases_df.file_name.str.contains('BRGeneva')
# plot_this(bases_df, idcs_lib, wave_bases, flux_bases)
#
# # Plot by library
# idcs_lib = (bases_df.z_star == 0.0190) & (bases_df.age_yr < 1.5e7)
# plot_this(bases_df, idcs_lib, wave_bases, flux_bases)
#
# # Plot by library
# idcs_lib = (bases_df.z_star == 0.0190) & (bases_df.age_yr > 5.5e6) & (bases_df.age_yr < 1.5e7)
# plot_this(bases_df, idcs_lib, wave_bases, flux_bases)




# # ------ Plot stellar bases age and metallicity
# fig, ax = plt.subplots(figsize=(9, 9))
#
# idcs_lib = (bases_df.age_yr < 1.5e7) & (bases_df.z_star > 0.005)
# z_lib_i = bases_df.loc[idcs_lib, 'z_star'].values
# age_lib_i = bases_df.loc[idcs_lib, 'age_yr'].values
# ax.scatter(np.log10(age_lib_i), z_lib_i, label='Group 1')
#
# ax.set_xlabel(r'$log(age)$')
# ax.set_ylabel(r'$z_{\star}$')
# ax.set_title('SSP bases metallicity versus age')
# ax.legend(loc='upper center')
# plt.show()
#



# # ------ Plot stellar bases age and metallicity
# fig, ax = plt.subplots(figsize=(9, 9))
# bases_root_list = ['Iku1', 'br04', 'BRGeneva', 'Mun1']
# for bases_lib in bases_root_list:
#
#     idcs_lib = bases_df['file_name'].str.contains(bases_lib)
#     if idcs_lib.any():
#         z_lib_i = bases_df.loc[idcs_lib, 'z_star'].values
#         age_lib_i = bases_df.loc[idcs_lib, 'age_yr'].values
#         ax.scatter(np.log10(age_lib_i), z_lib_i, label=bases_lib)
#
# ax.set_xlabel(r'$log(age)$')
# ax.set_ylabel(r'$z_{\star}$')
# ax.set_title('SSP bases metallicity versus age')
# ax.legend(loc='upper center')
# plt.show()


# # ------ Plot stellar bases age and metallicity
# SSP complete libraries
# param_min = 0.0
# for param_scale in ['Mcor_j', 'x_j']:
#
#     fig, ax = plt.subplots(figsize=(9, 9))
#     bases_root_list = ['Iku1', 'br04', 'BRGeneva', 'Mun1']
#     symbol_list = ['v', 'D', 'P', 'X']
#
#     print(param_scale)
#     for i_base, bases_lib in enumerate(bases_root_list):
#         idcs_lib = bases_df['file_name'].str.contains(bases_lib)
#         if idcs_lib.any():
#             z_lib_i = bases_df.loc[idcs_lib, 'z_star'].values
#             age_lib_i = bases_df.loc[idcs_lib, 'age_yr'].values
#             ax.scatter(np.log10(age_lib_i), z_lib_i, marker=symbol_list[i_base], label=bases_lib, color='tab:grey')
#
#     param_array = np.array(fit_output[param_scale])
#     age_array, z_array = np.array(fit_output['age_j']), np.array(fit_output['Z_j'])
#
#     idcs_ssp = param_array > param_min
#     scat = ax.scatter(np.log10(age_array[idcs_ssp]), z_array[idcs_ssp], c=param_array[idcs_ssp], cmap="RdYlGn", s=50, label='Measured populations')
#
#     ax.set_xlabel(r'$log(age)$')
#     ax.set_ylabel(r'$z_{\star}$')
#     ax.set_title('SSP bases metallicity versus age')
#     ax.legend(loc='upper center')
#
#     cbar = fig.colorbar(scat, orientation='horizontal')
#     cbar.set_label(param_scale)
#     cbar.formatter.useOffset = False
#     cbar.update_ticks()
#     plt.tight_layout()
#     plt.show()
