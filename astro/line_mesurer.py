import time
import numpy as np
import matplotlib.pyplot as plt
from tools.line_tools import LineMeasurer, gaussFunc

# from .measuring_tools.line_tools import gaussFunc, LineMeasurer
from scipy import stats, optimize, integrate
# Fake data

# Function selector
lm = LineMeasurer()

wave = np.linspace(4950, 5050)
m, n, noise = 0.0, 2.0, np.random.normal(0, 0.05, wave.size)
flux_cont = (m * wave + n) + noise
ampTrue, muTrue, sigmaTrue = 10, 5007, 2.3
flux_gauss = gaussFunc((wave, flux_cont), ampTrue, muTrue, sigmaTrue)
wave_regions = np.array([4960, 4980, 4996, 5015, 5030, 5045])
areaTrue = np.sqrt(2 * np.pi * sigmaTrue ** 2) * ampTrue

# Load observed spectrum
lm.load_spectrum(wave, flux_gauss)

# Declare regions data
idcsLines, idcsContinua = lm.define_masks(wave_regions)

# Identify line regions
lm.line_properties(idcsLines, idcsContinua, bootstrap_size=1000)

# Fit gaussian profile
lm.gauss_mcfit(idcsLines, idcsContinua, bootstrap_size=1000)

# Comparing flux integration techniques
lineWave, lineFlux = lm.wave[idcsLines], lm.flux[idcsLines]
continuaWave, continuaFlux = lm.wave[idcsContinua], lm.flux[idcsContinua]
lineContinuumFit = lineWave * lm.m_continuum + lm.n_continuum
areaSimps = integrate.simps(lineFlux, lineWave) - integrate.simps(lineContinuumFit, lineWave)
areaTrapz = integrate.trapz(lineFlux, lineWave) - integrate.trapz(lineContinuumFit, lineWave)
areaIntgPixel = (lm.flux[idcsLines].sum() - lineContinuumFit.sum()) * lm.pixelWidth
resampleWaveLine = np.linspace(lineWave[0]-10, lineWave[-1]+10, 100)
resampleFluxCont = resampleWaveLine * lm.m_continuum + lm.n_continuum
gaussianCurve = gaussFunc((resampleWaveLine, resampleFluxCont), *lm.p1)
areaGauss = (gaussianCurve.sum() - resampleFluxCont.sum()) * np.diff(resampleWaveLine).mean()

# Print the results
print(f'True area : {areaTrue}')
print(f'Simpsons rule: {areaSimps}')
print(f'Trapezoid rule: {areaTrapz}')
print(f'Pixel intgr: {areaIntgPixel}')
print(f'Pixel intgr MC: {lm.lineIntgFlux} +/- {lm.lineIntgErr}')
print(f'Gauss intgr: {areaGauss}')
print(f'Gauss intgr MC: {lm.lineGaussFlux} +/- {lm.lineGaussErr}')
print(f'True amp = {ampTrue}, mu = {muTrue}, sigma = {sigmaTrue}')
print(f'Fit amp = {lm.p1[0]:2f} +/- {lm.p1_Err[0]:2f}, mu = {lm.p1[1]:2f} +/- {lm.p1_Err[1]:2f}, '
      f'sigma = {lm.p1[2]:2f} +/- {lm.p1_Err[2]:2f}')
print(f'Eqw MC {lm.eqw} +/- {lm.eqwErr}')



# Graphical results
# print("Eqw ", self.lineIntgFlux / lineContinuumFit.mean(), self.lineIntgErr / lineContinuumFit.mean())
# print("Eqw MC ", eqw, eqwErr)


# print("Line Gauss Area", lineAreaSingle)
# print("Line Gauss Area MC", self.lineGaussFlux, self.lineGaussErr)

# print("Line MC integration", self.lineIntgFlux, self.lineIntgErrFlux)

# # Plot the results
# resampleWaveLine = np.linspace(lineWave[0]-10, lineWave[-1]+10, 100)
# resampleWaveCont = resampleWaveLine * slope + intercept
# gaussianCurve = gaussFunc((resampleWaveLine, resampleWaveCont), *fitParams)
# sumGaussian = (np.sum(gaussianCurve) - np.sum(resampleWaveCont)) * np.diff(resampleWaveLine).mean()
#


# y_array = np.array([1,2,3,2,1])
#
fig, ax = plt.subplots()
ax.plot(wave, flux_gauss, label='Observed line')
ax.scatter(continuaWave, continuaFlux, label='Continuum regions')
ax.plot(lineWave, lineContinuumFit, label='Observed line', linestyle=':')
ax.plot(resampleWaveLine, gaussianCurve, label='Gaussian fit', linestyle=':')
ax.legend()
ax.update({'xlabel':'Flux', 'ylabel':'Wavelength', 'title':'Gaussian fitting'})
plt.show()


# Code for line spectra sheet
# Declare data location
# objName = '8'
# root_folder = 'E:\\Dropbox\\Astrophysics\\Data\\WHT_observations\\bayesianModel\\'
# objectFolder = '{}{}\\'.format(root_folder, objName)
# dataLocation = '{}{}\\{}_objParams.txt'.format(root_folder, objName, objName)
#
# # Load the spectrum
# obsData = specS.load_obsData(dataLocation, objName)
# wave, flux = obsData['obs_wavelength'], obsData['obs_flux']
#
# # Normalizing the spectrum # TODO in here we use the norm flux from the continuum
# normFluxCoeff = np.median(flux)
# fluxNorm = flux / normFluxCoeff
#
# # Load lines data # TODO in here we use the lines proposed by the user
# objLinesLogDf = pd.read_csv(obsData['obj_lines_file'], delim_whitespace=True, header=0, index_col=0)
# idcs_validLines = ~objLinesLogDf.index.isin(['O2_7319A_b', 'H1_6563A_w'])
# linesDf = objLinesLogDf.loc[idcs_validLines]
# linesDf = linesDf.loc[:'H1_9546A']
# lineLabels = linesDf.index.values
# areaWaveN_matrix = linesDf.loc[:, 'w1':'w6'].values  # TODO in here we use the lines proposed by the user
#
# # Get line and adjacent continua region
# ares_indcs = np.searchsorted(wave, areaWaveN_matrix)
# idcsLines = (wave[ares_indcs[:, 2]] <= wave[:, None]) & (wave[:, None] <= wave[ares_indcs[:, 3]])
# idcsContinua = ((wave[ares_indcs[:, 0]] <= wave[:, None]) & (wave[:, None] <= wave[ares_indcs[:, 1]])) | (
#             (wave[ares_indcs[:, 4]] <= wave[:, None]) & (wave[:, None] <= wave[ares_indcs[:, 5]]))
#
#
# # Lines mask
# lines_mask = generate_object_mask(linesDf, wave, lineLabels)
#
# n_randomPoints = 1000
# rangeFittings = np.arange(n_randomPoints)
# n_lines = ares_indcs.shape[0]
# recombLinesIdx = linesDf.index.str.contains('He1') + linesDf.index.str.contains('He2') + linesDf.index.str.contains('H1')
# removeContinuumCheck = True
#
# p1_matrix = np.empty((n_randomPoints, 3))
# sqrt2pi = np.sqrt(2 * np.pi)
#
# fig, ax = plt.subplots(1, 1, figsize=(10, 6))
# for i in np.arange(n_lines):
#
#     # Get line_i wave and continua
#     lineWave, lineFlux = wave[idcsLines[:, i]], fluxNorm[idcsLines[:, i]]
#     continuaWave, continuaFlux = wave[idcsContinua[:, i]], fluxNorm[idcsContinua[:, i]]
#
#     #lineFit = lm.lineFit(wave, fluxNorm, idcsLines[:, i], plot_data=True)
#
#     # Compute linear line continuum and get the standard deviation on the continuum
#     slope, intercept, r_value, p_value, std_err = stats.linregress(continuaWave, continuaFlux)
#     continuaFit = continuaWave * slope + intercept
#     std_continuum = np.std(continuaFlux - continuaFit)
#     lineContinuumFit = lineWave * slope + intercept
#     continuumInt = lineContinuumFit.sum()
#     centerWave = lineWave[np.argmax(lineFlux)]
#     centerContInt = centerWave * slope + intercept
#     lineWidth = (lineWave[-1] - lineWave[0]) / lineWave.size
#     lineWidth2 = (lineWave[-1] - lineWave[0]) / (lineWave.size - 1)
#
#     # Compute matrix with random noise from the continua standard deviation
#     normalNoise = np.random.normal(0.0, std_continuum, (n_randomPoints, lineWave.size))
#     line_iFluxMatrix = lineFlux + normalNoise
#
#     # Compute integrated flux
#     areasArray = line_iFluxMatrix.sum(axis=1) - continuumInt
#     integInt, integStd = areasArray.mean(), areasArray.std()
#
#     # Initial values for fit
#     p0 = (lineFlux.max(), lineWave.mean(), 1.0)
#
#     print lineWave
#     print lineFlux
#     print continuaWave
#     print continuaFlux
#     print lineContinuumFit
#     print line_iFluxMatrix[0]
#     print
#     print
#
#     # Perform fit in loop
#     for j in rangeFittings: # This one is not powerfull enought... add more points
#         p1_matrix[j], pcov = optimize.curve_fit(gaussFunc, (lineWave, lineContinuumFit), line_iFluxMatrix[j], p0=p0)
#
#     # Compute mean values and std from gaussian fit
#     A, mu, sigma = p1_matrix[:, 0].mean(), p1_matrix[:, 1].mean(), p1_matrix[:, 2].mean()
#     gIntArray = p1_matrix[:, 0] * p1_matrix[:, 2] * sqrt2pi
#     p1Mean, gInt = p1_matrix.mean(axis=0), gIntArray.mean()
#     p1Std, gIntStd = p1_matrix.std(axis=0), gIntArray.std()
#
#     resampleGaus = np.linspace(lineWave[0]-10, lineWave[-1]+10, 100)
#     resampleContinuum = resampleGaus * slope + intercept
#     gaussianCurve = gaussFunc((resampleGaus, resampleContinuum), *p1Mean)
#     gaussianFlux2 = gaussFunc((resampleGaus, np.zeros(resampleGaus.size)), *p1Mean)
#     AreaCurve = simps(lineFlux, lineWave) - simps(lineContinuumFit, lineWave)
#
#     lineFit = lm.lineFit(wave, fluxNorm, idcsLines[:, i], idcsContinua[:, i], plot_data=True)
#
#     print lineLabels[i], lineFlux.sum() - continuumInt, (fluxNorm * lines_mask[i]).sum() - continuumInt, integInt, gInt, AreaCurve
#
#     ax.plot(lineWave, lineFlux, color='black', label='lines')
#     ax.plot(continuaWave, continuaFlux, '-', color='tab:blue', label='continuum')
#     ax.plot(continuaWave, continuaFit, ':', color='tab:green', label='fit continuum')
#     ax.plot(lineWave, lineContinuumFit, ':', color='tab:green', label='fit continuum raro')
#     ax.plot(resampleGaus, gaussianCurve, ':', color='tab:purple', label='Fit Gaussian')
#     ax.plot(lineFit['resampleRegion'], lineFit['resampleCurve'], ':', color='tab:red', label='New fit Gaussian')
#
#     # ax.plot(resampleGaus, gaussianFlux2, color='tab:cyan', label='gaussian fit2')
#
#     ax.update({'xlabel': r'Wavelength $(\AA)$', 'ylabel': 'Flux (normalised)'})
#     ax.legend()
#     plt.show()