import numpy as np
import pandas as pd
import lime

from pathlib import Path
from astropy.io import fits
from arviz.stats.diagnostics import rhat
from lime.plots import STANDARD_PLOT
from matplotlib import pyplot as plt, rcParams, colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Load the fitting configuration
obsData = lime.load_cfg('../muse_CGCG007.ini')
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
conf_fit = lime.load_cfg(dataFolder/'HII_CHIM_TRY_grid.cfg')
plots_folder = Path(obsData['data_location']['plots_folder'])

# Load the photoionization grid
file_address = f'{dataFolder}/HII-CHI-mistry_1Myr_grid_O.txt'
grid_3D_DF = pd.read_csv(file_address, delim_whitespace=True, header=0)

# Lines for the fitting
obs_lines = np.array(['O3_4959A', 'O3_5007A',
                      'He1_5876A',
                      'N2_6548A', 'H1_6563A', 'N2_6584A',
                      'S2_6716A', 'S2_6731A',
                      'O2_7319A', 'O2_7330A',
                      'S3_9069A'])

# # Prepare grid interpolators
# model_variables = ['logOH', 'logU', 'logNO']
# gw = sr.GridWrapper()
# grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_3D_DF, axes_columns=model_variables)
# grid_interpolators = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')
#
# grid_lines = np.array(list(grid_dict.keys()))
# idcs_lines = np.isin(obs_lines, grid_lines)
# input_lines = obs_lines[idcs_lines]

# Grid to perform the fitting
# logOH_range = np.round(np.linspace(7.15, 9.0, 7), 3)
# logU_range = np.round(np.linspace(-3.90, -1.40, 7), 3)
# logNO_range = np.round(np.linspace(-1.90, -0.01, 7), 3)

logOH_range = np.array([7.15, 7.458, 7.767, 8.075, 8.383, 8.692, 9.])
logU_range = np.array([-3.9, -3.483, -3.067, -2.65, -2.233, -1.817, -1.4])
logNO_range = np.array([-1.9, -1.585, -1.27, -0.955, -0.64, -0.325, -0.01])

# Output file for the fitting
outputFits = resultsFolder/'grid_logOH_logU_logNO_advi.fits'

# ----------------------------- Recovering the grid fits ---------------------------------------------------------------

i_step, n_steps = 0, logNO_range.size * logOH_range.size * logNO_range.size

logOH_fit = np.empty(n_steps)
logNO_fit = np.empty(n_steps)
logU_fit = np.empty(n_steps)
r_hat_logOH = np.empty(n_steps)
r_hat_logU = np.empty(n_steps)
r_hat_logNO = np.empty(n_steps)
r_hat_mean = np.empty(n_steps)
diff_mean_percent = np.empty(n_steps)

lines_dict = {}
plot_lines = ['N2_6584A', 'O3_5007A', 'S2_6731A', 'He1_5876A']
for line in plot_lines:
    lines_dict[line] = np.empty(n_steps)

with fits.open(outputFits) as hdul:
    for i, logOH in enumerate(logOH_range):
        for j, logU in enumerate(logU_range):
            for k, logNO in enumerate(logNO_range):

                # True value coordinate for interpolation
                ext_name = f'{logOH * 1000:.0f}_{logU * -1000:.0f}_{logNO * -1000:.0f}'

                data = hdul[f'{ext_name}_traces'].data
                hdr_inputs = hdul[f'{ext_name}_inputs'].header

                logOH_trace = data['logOH']
                logU_trace = data['logU']
                logNO_trace = data['logNO']

                # logOH_fit[i_step] = hdul[f'{ext_name}_outputs'].header['logOH']
                # logU_fit[i_step] = hdul[f'{ext_name}_outputs'].header['logU']
                # logNO_fit[i_step] = hdul[f'{ext_name}_outputs'].header['logNO']

                logOH_fit[i_step] = logOH
                logU_fit[i_step] = logU
                logNO_fit[i_step] = logNO

                r_hat_logU[i_step] = rhat(logU_trace.reshape(10, 500))
                r_hat_logOH[i_step] = rhat(logOH_trace.reshape(10, 500))
                r_hat_logNO[i_step] = rhat(logU_trace.reshape(10, 500))
                r_hat_mean[i_step] = (r_hat_logU[i_step] + r_hat_logOH[i_step] + r_hat_logNO[i_step])/3

                for line in plot_lines:
                    lines_dict[line][i_step] = data[line].mean()/hdr_inputs[line] - 1

                OH_diff = (1 - logOH_trace.mean()/logOH)
                NO_diff = (1 - logNO_trace.mean()/logNO)
                logU_diff = (1 - logU_trace.mean()/logU)
                diff_mean_percent[i_step] = (np.abs(OH_diff) + np.abs(NO_diff) + np.abs(logU_diff))/3 * 100

                i_step += 1

# ----------------------------- R_hat plot --------------------------------------------------------------------------


STANDARD_PLOT['axes.labelsize'] = 20
STANDARD_PLOT['ytick.labelsize'] = 16
STANDARD_PLOT['xtick.labelsize'] = 16

# Plot combined mask
defaultConf = STANDARD_PLOT.copy()
rcParams.update(defaultConf)

fig = plt.figure(figsize=(10, 10), dpi=600)
ax = fig.add_subplot(projection='3d')
sc = ax.scatter(logOH_fit, logNO_fit, logU_fit, c=r_hat_mean)
# sc = ax.scatter(logOH_fit, logNO_fit, logU_fit, c=lines_dict['N2_6584A'])

cbar = fig.colorbar(sc, fraction=0.025)
cbar.set_label(r'$\hat{R}$', rotation=0, labelpad=30)

ax.xaxis.set_major_locator(plt.MaxNLocator(4))
ax.yaxis.set_major_locator(plt.MaxNLocator(3))
ax.zaxis.set_major_locator(plt.MaxNLocator(4))

ax.set_xlabel(r'$log(O/H)$')
ax.set_ylabel(r'$log(N/O)$')
ax.set_zlabel(r'$log(U)$')
ax.xaxis.labelpad=20
ax.yaxis.labelpad=20
ax.zaxis.labelpad=20
ax.elev = 2.5
ax.azim = 276.5
# plt.show()
plt.tight_layout()
plt.savefig(plots_folder/'r_hat_tests_grid.png', bbox_inches='tight')

# ----------------------------- Parameter deviation --------------------------------------------------------------------------

# Plot combined mask
defaultConf = STANDARD_PLOT.copy()
rcParams.update(defaultConf)

fig = plt.figure(figsize=(10, 10), dpi=600)
ax = fig.add_subplot(projection='3d')
sc = ax.scatter(logOH_fit, logNO_fit, logU_fit, c=diff_mean_percent, norm=colors.LogNorm())

cbar = fig.colorbar(sc, fraction=0.025)
cbar.set_label(r'$\%$ mean deviation', labelpad=30)

ax.xaxis.set_major_locator(plt.MaxNLocator(4))
ax.yaxis.set_major_locator(plt.MaxNLocator(3))
ax.zaxis.set_major_locator(plt.MaxNLocator(4))

ax.set_xlabel('log(O/H)')
ax.set_ylabel('log(N/O)')
ax.set_zlabel('logU')
ax.xaxis.labelpad=20
ax.yaxis.labelpad=20
ax.zaxis.labelpad=20
ax.elev = 2.5
ax.azim = 276.5
# plt.show()

plt.tight_layout()
plt.savefig(plots_folder/'percentage_discrepancy.png', bbox_inches='tight')