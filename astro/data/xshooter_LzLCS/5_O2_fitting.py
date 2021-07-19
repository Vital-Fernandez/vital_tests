import numpy as np
from astropy.io import fits
from pathlib import Path
import src.specsiser as sr
from matplotlib import pyplot as plt
import pymc3 as pm
import arviz as az
import theano.tensor as tt

def gaussian_model(x, amp, center, sigma):
    """1-d gaussian curve : gaussian(x, amp, cen, wid)"""
    return amp * np.exp(-0.5 * (((x-center)/sigma) * ((x-center)/sigma)))

def mixture_density_mult(w, mu, sd, x):
    logp = tt.log(w) + pm.Normal.dist(mu, sd).logp(x)
    return tt.sum(tt.exp(logp), axis=1)

obsData = sr.loadConfData('./xshooter_LzLCS.ini')
data_folder = Path(obsData['data_location']['data_folder'])
results_folder = Path(obsData['data_location']['results_folder'])
objfile_list = obsData['data_location']['objfile_list']
sigmafile_list = obsData['data_location']['sigmafile_list']
objRef_list = obsData['data_location']['ref_list']
maskfile = obsData['data_location']['generalMask']

wmin_array = obsData['sample_data']['w_min_array']
wmax_array = obsData['sample_data']['w_max_array']
norm_flux = obsData['sample_data']['norm_flux']
z_obj = obsData['sample_data']['z_obj']
profile_conf = obsData['line_fitting']

verbose = False

# Input data
i = 0
objName = objRef_list[i]
spec_file, sigm_file = data_folder/objfile_list[i], data_folder/sigmafile_list[i]

# Output data
lineslog_file = results_folder/f'{objName}_linesLog.txt'
lineslog_table = results_folder/f'{objName}_flux_table'

# Load inputs
wave, flux, header = sr.import_fits_data(spec_file, instrument='xshooter', frame_idx=0)
wave_sigma, sigma, header_sigma = sr.import_fits_data(sigm_file, instrument='xshooter', frame_idx=0)

lm = sr.LineMesurer(wave, flux, crop_waves=[wmin_array[i], wmax_array[i]], input_err=sigma, normFlux=norm_flux, redshift=z_obj)

mask_local = data_folder/f'{objName}_mask'
mask_local_df = sr.lineslogFile_to_DF(mask_local)

lineLabel = 'O2_3726A_b'
wave_regions = mask_local_df.loc[lineLabel, 'w1':'w6'].values
lm.fit_from_wavelengths(lineLabel, wave_regions, user_conf=profile_conf)
lm.print_results(show_plot=True, show_fit_report=True, log_scale=False)

# Input data
idcsEmis, idcsCont = lm.define_masks(lm.wave_rest, lm.flux, wave_regions)

idcsLine = idcsEmis + idcsCont
x_array = lm.wave[idcsLine]
y_array = lm.flux[idcsLine]
sigma_array = lm.errFlux[idcsLine]
w_array = 1.0 / lm.errFlux[idcsLine] if lm.errFlux is not None else np.full(x_array.size, 1.0 / lm.std_cont)
y_lineal = lm.m_cont * x_array + lm.n_cont

# plt.step(x_array, y_array)
# plt.plot(x_array, y_lineal)
# plt.plot(x_array, lm.errFlux[idcsLine])
# plt.legend()
# plt.show()

n_comps = 4
n_pixels = len(y_lineal)
range_comps = np.arange(n_comps)
mu_theo = np.array([3726.0, 3729.0, 3726.0, 3729.0]) * (1 + 0.282931)


x_in = x_array[:, None]
y_in = y_array[:, None]
#-----------------------------------------------------------------------------


# with pm.Model():
#
#     # Model Priors
#     amp_array = pm.HalfNormal('amp_array', 50., shape=n_comps)
#     my_array = pm.Normal('mu_array', mu_theo, 3., shape=n_comps)
#     sigma_array = pm.HalfCauchy('sigma_array', 5., shape=n_comps)
#     pixelNoise = pm.HalfCauchy('pixelNoise', 0.5)
#
#     # Theoretical line profiles
#     theoFlux_i = mixture_density_mult(amp_array, my_array, sigma_array, x_in) + y_lineal
#
#     # Model likelihood
#     pm.Normal('emission_Y', theoFlux_i, pixelNoise, observed=y_array)
#
#     # Run sampler
#     trace = pm.sample(draws=1000, tune=1000, chains=2, cores=2)
#
# pm.summary(trace)
#
# amp_trace, mu_trace, sigma_trace = trace['amp_array'], trace['mu_array'], trace['sigma_array']
# amp_mean, mu_mean, sigma_mean = amp_trace.mean(axis=0), mu_trace.mean(axis=0), sigma_trace.mean(axis=0)
# print(trace['amp_array'].mean(axis=0))
#
# wave_resample = np.linspace(x_array[0], x_array[-1], x_array.size * 20)
# cont_resample = lm.m_cont * wave_resample + lm.n_cont
# hmc_curve = mixture_density_mult(amp_mean, mu_mean, sigma_mean, wave_resample[:, None]).eval()
#
# # Plot the results
# fig, ax = plt.subplots()
# ax.step(x_array, y_array, label='Object spectrum')
# ax.plot(wave_resample, hmc_curve + cont_resample, label='Fitting result',  color='tab:red')
# for i in range(n_comps):
#     hmc_curve_i = mixture_density_mult(amp_mean[i], mu_mean[i], sigma_mean[i], wave_resample[:, None]).eval()
#     ax.plot(wave_resample, hmc_curve_i + cont_resample, label='Fitting result', linestyle=':')
# ax.legend()
# ax.update({'xlabel':r'Wavelength $(\AA)$', 'ylabel':r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1}) \cdot 10^{17}$',
#            'title':'Gaussian fitting'})
# plt.show()


#-----------------------------------------------------------------------------

with pm.Model() as model3:
    amp = pm.Uniform('amp', 0, 100, shape=n_comps)
    loc = pm.Normal('mu', mu_theo, 3, shape=n_comps)
    scale = pm.Uniform('scale', 0, 10, shape=n_comps)

    flux_comp = np.zeros(n_pixels)
    flux_comp = flux_comp + y_lineal
    for i in range_comps:
        flux_comp = flux_comp + gaussian_model(x_array, amp[i], loc[i], scale[i])
    # gauss = pm.Deterministic('gauss', flux_comp)

    y = pm.Normal('y', mu=flux_comp, sigma=sigma_array, observed=y_array)

    trace = pm.sample(2000, tune=2000, chains=3, cores=3, init='advi')


amp_trace, mu_trace, sigma_trace = trace['amp'], trace['mu'], trace['scale']
amp_mean, mu_mean, sigma_mean = amp_trace.mean(axis=0), mu_trace.mean(axis=0), sigma_trace.mean(axis=0)
print(trace['amp'].mean(axis=0))

wave_resample = np.linspace(x_array[0], x_array[-1], x_array.size * 20)
cont_resample = lm.m_cont * wave_resample + lm.n_cont
hmc_curve_total = np.zeros(wave_resample.size)
for i in range(n_comps):
    hmc_curve_total += gaussian_model(wave_resample, amp_mean[i], mu_mean[i], sigma_mean[i])

# Plot the results
fig, ax = plt.subplots()
ax.step(x_array, y_array, label='Object spectrum')
ax.plot(wave_resample, hmc_curve_total + cont_resample, label='Fitting result',  color='tab:red')
for i in range(n_comps):
    hmc_curve_i =  gaussian_model(wave_resample, amp_mean[i], mu_mean[i], sigma_mean[i])
    ax.plot(wave_resample, hmc_curve_i + cont_resample, label='Fitting result', linestyle=':')
ax.legend()
ax.update({'xlabel':r'Wavelength $(\AA)$', 'ylabel':r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1}) \cdot 10^{17}$',
           'title':'Gaussian fitting'})
plt.show()

plt.step(x_array, y_array)
plt.plot(x_array, y_lineal)
plt.plot(x_array, lm.errFlux[idcsLine])
plt.legend()
plt.show()




# import numpy as np
# from pymc3 import *
# import matplotlib.pyplot as plt
#
# # set random seed for reproducibility
# np.random.seed(12345)
#
# x = np.arange(5,400,10)*1e3
#
# # Parameters for gaussian
# amp_true = 0.2
# size_true = 1.8
# ps_true = 0.1
#
# #Gaussian function
# gauss = lambda x,amp,size,ps: amp*np.exp(-1*(np.pi**2/(3600.*180.)*size*x)**2/(4.*np.log(2.)))+ps
# f_true = gauss(x=x,amp=amp_true, size=size_true, ps=ps_true )
#
# # add noise to the data points
# noise = np.random.normal(size=len(x)) * .02
# f = f_true + noise
# f_error = np.ones_like(f_true)*0.05*f.max()
#
# with Model() as model3:
#     amp = Uniform('amp', 0.05, 0.4, testval= 0.15)
#     size = Uniform('size', 0.5, 2.5, testval= 1.0)
#     ps = Normal('ps', 0.13, 40, testval=0.15)
#
#     gauss = Deterministic('gauss', amp*np.exp(-1*(np.pi**2*size*x/(3600.*180.))**2/(4.*np.log(2.)))+ps)
#
#     y =Normal('y', mu=gauss, tau=1.0/f_error**2, observed=f)
#
#     start=find_MAP()
#     step=NUTS()
#     trace=sample(2000, start=start)
#
# # extract and plot results
# y_min = np.percentile(trace.gauss,2.5,axis=0)
# y_max = np.percentile(trace.gauss,97.5,axis=0)
# y_fit = np.percentile(trace.gauss,50,axis=0)
# plt.plot(x, f_true, 'b', marker='None', ls='-', lw=1, label='True')
# plt.errorbar(x, f, yerr=f_error, color='r', marker='.', ls='None', label='Observed')
# plt.plot(x, y_fit, 'k', marker='+', ls='None', ms=5, mew=1, label='Fit')
# plt.fill_between(x, y_min, y_max, color='0.5', alpha=0.5)
# plt.legend()
# plt.show()


