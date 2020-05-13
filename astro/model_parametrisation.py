# import numpy as np
# from src.specsyzer.physical_model.line_tools import LineMeasurer, gaussFunc
# from scipy import integrate
# from matplotlib import pyplot as plt, rcParams
# import lmfit
#
# lm = LineMeasurer()
#
# # Generate fake data
# wave = np.linspace(4950, 5050, 300)
# m, n, noise = 0.0, 2.0, np.random.normal(0, 0.05, wave.size)
# zero_level = (m * wave + n) #+ noise
# ampTrue, muTrue, sigmaTrue = 10, 5007, 2.3
# fluxObs = gaussFunc((wave, zero_level), ampTrue, muTrue, sigmaTrue)
# wave_regions = np.array([4960, 4980, 4996, 5015, 5030, 5045])
# areaTrue = np.sqrt(2 * np.pi * sigmaTrue ** 2) * ampTrue
#
# fluxObsLog = np.log(fluxObs)
# lm = LineMeasurer(wave, fluxObsLog)
#
# # Declare regions data
# idcsLines, idcsContinua = lm.define_masks(wave_regions)
#
# # Identify line regions
# lm.line_properties(idcsLines, idcsContinua, bootstrap_size=1000)
#
# # Fit gaussian profile
# lm.line_fitting(idcsLines, idcsContinua, bootstrap_size=1000)
#
# # Comparing flux integration techniques
# lineWave, lineFlux = lm.wave[idcsLines], lm.flux[idcsLines]
# continuaWave, continuaFlux = lm.wave[idcsContinua], lm.flux[idcsContinua]
# resampleWaveLine = np.linspace(lineWave[0] - 10, lineWave[-1] + 10, 100)
# resampleFluxCont = resampleWaveLine * lm.m_continuum + lm.n_continuum
# gaussianCurve = gaussFunc((resampleWaveLine, resampleFluxCont), *lm.p1)
# linearFitCont = lm.wave * lm.m_continuum + lm.n_continuum
# logGaussCurve = (np.log(ampTrue) + (-0.5 * np.power((lm.wave-muTrue)/sigmaTrue, 2))) + linearFitCont
#
#
# # -------------- LMFIT --------------
# from lmfit import Model
#
# def gaussian(x, amp, cen, wid):
#     """1-d gaussian: gaussian(x, amp, cen, wid)"""
#     return (amp / (np.sqrt(2*np.pi) * wid)) * np.exp(-(x-cen)**2 / (2*wid**2))
#
#
# def line(x, slope, intercept):
#     """a line"""
#     return slope*x + intercept
#
# mod = Model(gaussian) + Model(line)
# pars = mod.make_params(amp=lm.flux.max(), cen=lm.wave[np.argmax(lm.flux)], wid=2,
#                        slope=0, intercept=1)
# result = mod.fit(lm.flux, pars, x=lm.wave)
# print(result.fit_report())
#
# # -------------- PLOT RESULTS --------------
# fig, ax = plt.subplots()
# ax.scatter(lm.wave, lm.flux, label='Observed line')
# ax.plot(resampleWaveLine, gaussianCurve, label='Gaussian fit', linestyle=':')
# ax.plot(lm.wave, linearFitCont, label='Linear fit')
# ax.plot(lm.wave, result.best_fit, label='LMFIT')
# ax.plot(lm.wave, logGaussCurve, label='Log Gauss Curve')
#
# ax.legend()
# ax.update({'xlabel': 'Flux', 'ylabel': 'Wavelength', 'title': 'Gaussian fitting'})
# plt.show()



# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.optimize import curve_fit
#
#
# def gauss(x, H, A, x0, sigma):
#     return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))
#
# def gauss_fit(x, y):
#     mean = sum(x * y) / sum(y)
#     sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
#     popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
#     return popt
#
#
# # generate simulated data
# np.random.seed(123)  # comment out if you want different data each time
# xdata = np.linspace(3, 10, 100)
# ydata_perfect = np.log10(gauss(xdata, 20, 5, 6, 1))
# ydata = ydata_perfect
#
# H, A, x0, sigma = gauss_fit(xdata, ydata)
# FWHM = 2.35482 * sigma
#
# print('The offset of the gaussian baseline is', H)
# print('The center of the gaussian fit is', x0)
# print('The sigma of the gaussian fit is', sigma)
# print('The maximum intensity of the gaussian fit is', H + A)
# print('The Amplitude of the gaussian fit is', A)
# print('The FWHM of the gaussian fit is', FWHM)
#
# plt.plot(xdata, ydata, 'ko', label='data')
# plt.plot(xdata, ydata_perfect, '-k', label='data (without_noise)')
# plt.plot(xdata, gauss(xdata, *gauss_fit(xdata, ydata)), '--r', label='fit')
#
#
# plt.legend()
# plt.title('Gaussian fit,  $f(x) = A e^{(-(x-x_0)^2/(2sigma^2))}$')
# plt.xlabel('Motor position')
# plt.ylabel('Intensity (A)')
# plt.show()
#
# print('True values', 20, 5, 6, 1)
#
# '''Using scipy.curve_fit to fit a normal distribution to data.'''

# import random
#
# import numpy as np
# from scipy.optimize import curve_fit
# import matplotlib.pyplot as plt
#
#
# # Create a function which returns a Gaussian (normal) distribution.
# def gauss(x, *p):
#     a, b, c, d = p
#     y = a*np.exp(-np.power((x - b), 2.)/(2. * c**2.)) + d
#
#     return y
#
# # Choose some starting params for our distribution and perturb them
# # using random numbers.
# p_initial = [1.0, 0.0, 0.1, 0.0]
# p_perturbed = p_initial
#
# N = 100 # Number of data points.
#
# # Create our data sets. Perturb the y-data with randomness and
# # generate completely random data for the errors.
# x = np.linspace(-1, 1, N)
# y = np.log10(gauss(x, *p_perturbed))
# e = np.array([random.random()*0.1 for _ in y])
#
#
#
# # Use curve_fit to fit the gauss function to our data. Use the
# # unperturbed p_initial as our initial guess.
# popt, pcov = curve_fit(gauss, x, y, p0=p_initial)
#
# # Generate y-data based on the fit.
# y_fit = gauss(x, *popt)
#
# # Create a plot of our work, showing both the data and the fit.
# fig, ax = plt.subplots()
#
# ax.errorbar(x,y,e)
# ax.plot(x, y_fit, color = 'red')
#
# ax.set_xlabel(r'$x$')
# ax.set_ylabel(r'$f(x)$')
# ax.set_title('Using scipy.curve_fit to fit a normal distribution.')
#
# plt.show()

import numpy as np
from src.specsyzer.physical_model.line_tools import LineMeasurer, gaussFunc
from scipy import integrate
from matplotlib import pyplot as plt, rcParams
import lmfit

lm = LineMeasurer()

# Generate fake data
wave = np.linspace(4950, 5050, 300)
m, n, noise = 0.0, 2.0, np.random.normal(0, 0.05, wave.size)
zero_level = (m * wave + n) #+ noise
ampTrue, muTrue, sigmaTrue = 10, 5007, 2.3
fluxObs = gaussFunc((wave, zero_level), ampTrue, muTrue, sigmaTrue)
wave_regions = np.array([4960, 4980, 4996, 5015, 5030, 5045])
areaTrue = np.sqrt(2 * np.pi * sigmaTrue ** 2) * ampTrue

equivFLux = np.log(ampTrue) - 0.5 * np.power((wave - muTrue)/sigmaTrue, 2)
fluxLine = gaussFunc((wave, np.zeros(wave.size)), ampTrue, muTrue, sigmaTrue)

flux = 0.850 #np.sqrt(2 * np.pi * sigmaTrue ** 2) * ampTrue
flux_err = flux * 0.05

logFlux = np.log10(flux)
logFluxErr = np.log10(flux + flux_err) - np.log10(flux)

fig, ax = plt.subplots()
ax.plot(wave, fluxObs, label='Gauss Curve')
ax.plot(wave, zero_level, ':', label='Zero level')
ax.plot(wave, np.log(fluxObs), label='log Gauss Curve')
ax.plot(wave, np.log(zero_level), ':', label='log Zero level')

# ax.plot(wave, fluxLine, label='Flux line')
# ax.plot(wave, np.log10(fluxLine), label='Flux line')
# ax.plot(wave, equivFLux, label='Flux Equiv')

ax.legend()
ax.update({'xlabel': 'Flux', 'ylabel': 'Wavelength', 'title': 'Gaussian fitting'})
plt.show()