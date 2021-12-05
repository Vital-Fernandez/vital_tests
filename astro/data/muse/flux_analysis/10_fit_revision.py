import numpy as np
import pandas as pd
from pathlib import Path
import src.specsiser as sr
from src.specsiser.inference_model import fits_db
from astro.data.muse.common_methods import grid_columns
from timeit import default_timer as timer
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams
from lmfit import models
from src.specsiser.data_printing import DARK_PLOT, background_color, foreground_color

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

# Load the data
grid_file = Path('D:/Dropbox/Astrophysics//Data/muse-Ricardo/Data/HII-CHI-mistry_1Myr_grid.csv')
grid_DF = pd.read_csv(grid_file, skiprows=1, names=grid_columns.values())
grid_DF.logNO = np.round(grid_DF.logNO.values, decimals=3)
df_columns = np.array(['wavelength', 'ion',	'intg_flux', 'intg_err'])

# Remove carbon dimension
idcs_rows = grid_DF.carbon == 'O'
idcs_columns = ~grid_DF.columns.isin(['carbon'])
grid_3D_DF = grid_DF.loc[idcs_rows, idcs_columns]

model_variables = ['logOH', 'logU', 'logNO']
gw = sr.ModelGridWrapper()
grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_3D_DF, axes_columns=model_variables)
grid_interpolators = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')

exclude_lines = np.array(['H1_4861A', 'N1_5198A', 'N1_5200A', 'C2_4267A', 'O2_4651A', 'C2_4659A'])
min_wavelength = 4700
input_lines = np.array(list(grid_dict.keys()))

ion_array, wave_array, latex_array = sr.label_decomposition(input_lines)
idcs_lines = ~np.isin(input_lines, exclude_lines) & (min_wavelength < wave_array)
input_lines = input_lines[idcs_lines]

r_steps = 5
logOH_range = np.round(np.linspace(7.15, 9.0, r_steps), 3)
logU_range = np.round(np.linspace(-3.90, -1.40, r_steps), 3)
logNO_range = np.round(np.linspace(-1.90, -0.01, r_steps), 3)
i_step, n_steps = 0, logNO_range.size * logOH_range.size * logNO_range.size

# Loop throught the grid of synthetic conditions (around 30 seconds per fit)
objFolder = resultsFolder / f'CGCG007/'
# outputFits = objFolder / f'grid_logOH_logU_logNO.fits'
outputFits = objFolder / f'grid_logOH_logU_logNO_advi.fits'

image_logOH = np.full((r_steps, r_steps), np.nan)
image_logU = np.full((r_steps, r_steps), np.nan)

t_steps = r_steps * r_steps * r_steps
logOH_fit = np.empty(t_steps)
logNO_fit = np.empty(t_steps)
logU_fit = np.empty(t_steps)
r_hat_logOH = np.empty(t_steps)
r_hat_logU = np.empty(t_steps)
r_hat_logNO = np.empty(t_steps)
r_hat_mean = np.empty(t_steps)
cHbeta_true = 0.255
err_mean = np.empty(t_steps)

# with fits.open(outputFits) as hdul:
#     for i, logOH in enumerate(logOH_range):
#         for j, logU in enumerate(logU_range):
#             for k, logNO in enumerate(logNO_range):
#
#                 print(f'- Fit {i_step}/{n_steps}: expected time {n_steps*30/3600:0.3f} hours')
#
#                 # True value coordinate for interpolation
#                 ext_name = f'{logOH*1000:.0f}_{logU*-1000:.0f}_{logNO*-1000:.0f}'
#
#                 logOH_measured[i_step] = hdul[f'{ext_name}_outputs'].header['logOH']
#                 logU_measured[i_step] = hdul[f'{ext_name}_outputs'].header['logU']
#                 logNO_measured[i_step] = hdul[f'{ext_name}_outputs'].header['logNO']
#
#                 logOH_fit[i_step] = hdul[f'{ext_name}_outputs'].header['logOH']
#                 logU_fit[i_step] = hdul[f'{ext_name}_outputs'].header['logU']
#                 logNO_fit[i_step] = hdul[f'{ext_name}_outputs'].header['logNO']
#                 cHbeta_fit[i_step] = hdul[f'{ext_name}_outputs'].header['cHbeta']
#
#                 err_array = np.array([1 - logOH_fit/logOH,
#                                       1 - logU_fit / logOH,
#                                       1 - logNO_fit / logOH,
#                                       1 - cHbeta_fit / cHbeta_true])
#
#
#                 r_hat_logU[i_step] = hdul[f'{ext_name}_traces'].header['logU_r_hat']
#                 r_hat_logOH[i_step] = hdul[f'{ext_name}_traces'].header['logOH_r_hat']
#                 r_hat_logNO[i_step] = hdul[f'{ext_name}_traces'].header['logOH_r_hat']
#                 r_hat_mean[i_step] = (r_hat_logU[i_step] + r_hat_logOH[i_step] + r_hat_logNO[i_step]) / 3.0
#                 err_mean[i_step] = err_array.mean()
#                 i_step += 1
#
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot(projection='3d')
# sc = ax.scatter(logOH_fit, logU_fit, logNO_fit, c=err_mean)
# cbar = fig.colorbar(sc)
# ax.set_xlabel('logOH')
# ax.set_ylabel('logU')
# ax.set_zlabel('logNO')
# plt.show()


with fits.open(outputFits) as hdul:
    for i, logOH in enumerate(logOH_range):
        for j, logU in enumerate(logU_range):
            for k, logNO in enumerate(logNO_range):

                print(f'- Fit {i_step}/{n_steps}: expected time {n_steps*30/3600:0.3f} hours')

                # True value coordinate for interpolation
                ext_name = f'{logOH*1000:.0f}_{logU*-1000:.0f}_{logNO*-1000:.0f}'

                # logOH_fit[i_step] = hdul[f'{ext_name}_outputs'].header['logOH']
                # logU_fit[i_step] = hdul[f'{ext_name}_outputs'].header['logU']
                # logNO_fit[i_step] = hdul[f'{ext_name}_outputs'].header['logNO']
                logOH_fit[i_step] = logOH
                logU_fit[i_step] = logU
                logNO_fit[i_step] = logNO

                r_hat_logU[i_step] = hdul[f'{ext_name}_traces'].header['logU_r_hat']
                r_hat_logOH[i_step] = hdul[f'{ext_name}_traces'].header['logOH_r_hat']
                r_hat_logNO[i_step] = hdul[f'{ext_name}_traces'].header['logOH_r_hat']
                r_hat_mean[i_step] = (r_hat_logU[i_step] + r_hat_logOH[i_step] + r_hat_logNO[i_step]) / 3.0
                i_step += 1

# Plot combined mask
defaultConf = DARK_PLOT.copy()
rcParams.update(defaultConf)

fig = plt.figure(figsize=(10, 10))
fig.patch.set_facecolor(background_color)
ax = fig.add_subplot(projection='3d')
sc = ax.scatter(logOH_fit, logU_fit, logNO_fit, c=r_hat_logOH)
cbar = fig.colorbar(sc,fraction=0.046)
cbar.set_label(r'R-hat', rotation=270, labelpad=30)

ax.set_xlabel('log(O/H)')
ax.set_ylabel('logU')
ax.set_zlabel('log(N/O)')
ax.xaxis.labelpad=15
ax.yaxis.labelpad=15
ax.zaxis.labelpad=15
plt.show()




