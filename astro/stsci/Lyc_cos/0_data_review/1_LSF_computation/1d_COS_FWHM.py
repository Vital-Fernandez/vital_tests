from pathlib import Path
import lime
from matplotlib import pyplot as plt, rc_context
import numpy as np
import re
from astro.stsci.tools import gaussian, fit_func
from scipy import optimize
from astropy import modeling
from astropy.io import ascii



# Data location
project_folder = Path('/')
obs_folder = Path('/home/vital/Astrodata/STScI')
lsf_folder =  Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_file')
lsf_fittings = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_fitted')

# Cfg file
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v0.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filecode"].str.contains(pattern, na=False)]

# Minimizers
fitter = modeling.fitting.LevMarLSQFitter()
model = modeling.models.Gaussian1D(amplitude=0.08, stddev=10.)

# Containers
gauss_LM = []
gauss_scipy = []
poly_scipy = []

# Loop thought the files
lsf_files = list(lsf_folder.glob('*.dat'))
for j, pname in enumerate(lsf_files):
    if j >= 0:
        print(j, pname)
        matrix_arr = np.loadtxt(pname)
        wave_arr = matrix_arr[0,:]
        lsf_arr = matrix_arr[1:, :]

        # Containers for the data
        gauss_LM = np.empty(lsf_arr.shape[1])
        gauss_scipy = np.empty(lsf_arr.shape[1])
        poly_scipy = np.empty(lsf_arr.shape[1])

        for i in np.arange(lsf_arr.shape[1]):
            x_arr = np.arange(lsf_arr.shape[0])
            fitted_model = fitter(model, x_arr, lsf_arr[:, i])

            ### scipy ###
            popt, _ = optimize.curve_fit(gaussian, x_arr, lsf_arr[:, i], p0=[0.05, x_arr.size/2, 4.])

            ### Gaussfit in python ###
            pars, cov = optimize.curve_fit(fit_func, x_arr, lsf_arr[:, i],
                                           p0=[0.1, x_arr.size / 2, 8., 0.00000481032, 8.11350e-07, -2.41179e-06],
                                           maxfev=100000)

            gauss_LM[i] = 2.35482 * np.absolute(fitted_model.stddev)
            gauss_scipy[i] = 2.35482 * np.absolute(popt[2])
            poly_scipy[i] = 2.35482 * pars[2]

            # # Compare the fittings
            # with rc_context({'figure.dpi': 250, 'figure.figsize': (4, 4)}):
            #     fig, ax = plt.subplots()
            #     ax.scatter(x_arr, lsf_arr[:, i], label=f'Wavelength {wave_arr[i]}Å LSF')
            #     ax.plot(x_arr, fitted_model(x_arr), '-', label=f'Levenberg-Marquardt Gaussian {gauss_LM[i]:.3f}', color='tab:red')
            #     ax.plot(x_arr, gaussian(x_arr, *popt), '--', label=f'Scipy Gaussian {gauss_scipy[i]:.3f}', color='tab:green')
            #     ax.plot(x_arr, fit_func(x_arr, *pars), ':', label=f'Scipy polynomial {poly_scipy[i]:.3f}', color='tab:orange')
            #     ax.legend(fontsize=6)
            #     plt.show()

        # # Show the FWHM evolution with wavelength
        # print(f'FWHM trend: {pname.name} ')
        # with rc_context({'figure.dpi': 250, 'figure.figsize': (4, 4)}):
        #     fig, ax = plt.subplots()
        #     ax.plot(wave_arr, gauss_LM, label=f'Levenberg-Marquardt Gaussian')
        #     ax.plot(wave_arr, gauss_scipy, label=f'Scipy Gaussian', linestyle='--')
        #     ax.plot(wave_arr, poly_scipy, label=f'Scipy polynomial', linestyle='--')
        #     ax.update({'xlabel': f'Wavelength (Å)', 'ylabel': f'FWHM (pixels)', 'title': pname.stem})
        #     ax.legend()
        #     plt.show()

        # Compute the polynomial
        fit  = np.poly1d(np.polyfit(wave_arr, np.array(gauss_scipy), 4))
        wave_all = wave_arr
        fwhm_arr_fit = fit(wave_all)

        with rc_context({'figure.dpi': 250, 'figure.figsize': (4, 4)}):
            fig, ax = plt.subplots()
            ax.plot(wave_arr, poly_scipy, label='Polynomical')
            ax.plot(wave_all, gauss_LM, label='gauss_LM')
            ax.plot(wave_all, fwhm_arr_fit, 'r', label='Interpolation', linestyle='--')
            ax.update({'xlabel': f'Wavelength (Å)', 'ylabel': f'FWHM (pixels)'})
            ax.legend()
            plt.show()

        ### save to ascii file ####
        ascii.write([wave_all, fwhm_arr_fit], lsf_fittings/f'{pname.stem}_fitted.dat', format='no_header',
                    overwrite=True)