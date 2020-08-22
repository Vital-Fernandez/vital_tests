import lmfit
import numpy as np
import pandas as pd
import pyneb as pn
import src.specsiser as sr
import matplotlib.pyplot as plt
from pathlib import Path
from src.specsiser.physical_model.gasContinuum_functions import NebularContinua


# Load data
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
files_path = Path('D:/Google drive/Astrophysics/Datos/MUSE - Amorin')
spectrum_file = 'voxel_cluster.txt'
idx_voxel = (170, 170)
flux_norm = 1e20
wave, flux_voxel = np.loadtxt(files_path/spectrum_file, unpack=True)
flux_voxel = flux_voxel / flux_norm

# Treatment tools
nebCalc = NebularContinua()
lm = sr.LineMeasurer(wave, flux_voxel)

# Remove the continuum
flux_noContinuum = lm.continuum_remover(noiseRegionLims=(5600, 5850))

noiseWaveLim = (5600, 5850)
intEmThreshold, intAbsThreshold = 1.5, 4

idcs_noiseRegion = (noiseWaveLim[0] <= wave) & (wave <= noiseWaveLim[1])
noise_mean, noise_std = flux_voxel[idcs_noiseRegion].mean(), flux_voxel[idcs_noiseRegion].std()
highLimit = intEmThreshold * (noise_mean + noise_std)
lowLimit = (noise_mean + noise_std) / intAbsThreshold

# Identify high flux regions # TODO Expand limit to absorption lines

idcsLineMask = np.where((flux_voxel >= lowLimit) & (flux_voxel <= highLimit))
wave_masked, flux_masked = wave[idcsLineMask], flux_voxel[idcsLineMask]

# Measure halpha flux
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
listLabels = ['H1_4861A', 'H1_6563A']
flux_dict = {}
for lineLabel in listLabels:
    wave_regions = linesDb.loc[lineLabel, 'w1':'w6'].values
    idcsLinePeak, idcsContinua = lm.define_masks(wave_regions)
    lm.line_properties(idcsLinePeak, idcsContinua, bootstrap_size=500)
    lm.gaussian_mcfit(idcsLinePeak, idcsContinua, bootstrap_size=500)
    flux_dict[lineLabel] = lm.lineIntgFlux

# Reddening correction
rc = pn.RedCorr()
rc.law = 'G03 LMC'
halpha_norm = flux_dict['H1_6563A']/flux_dict['H1_4861A']
rc.setCorr(obs_over_theo=halpha_norm, wave1=6563., wave2=4861.)
corr_spec = rc.getCorr(wave)
corr_Halpha = rc.getCorr(6563)
voxel_int = flux_voxel * corr_spec * flux_norm
voxel_int_masked = voxel_int[idcsLineMask]
halpha_flux = flux_dict['H1_6563A'] * corr_Halpha * flux_norm

# Compute nebular continuum
Te, ne = 10000.0, 100.0
HeII_HII, HeIII_HII = 0.1, 0.001
neb_int = nebCalc.flux_spectrum(wave, Te, halpha_flux, HeII_HII, HeIII_HII)

# Fit the nebular component
nebModel = lmfit.Model(nebCalc.flux_spectrum, prefix='neb_')
nebModel.set_param_hint('Halpha_Flux', value=halpha_flux, min=halpha_flux*0.90, max=halpha_flux*1.10)
nebModel.set_param_hint('Te', value=Te, min=8000, max=20000, vary=False)
nebModel.set_param_hint('He1_abund', value=HeII_HII, vary=False)
nebModel.set_param_hint('He2_abund', value=HeIII_HII, vary=False)
nebParams = nebModel.make_params()
out_neb = nebModel.fit(voxel_int_masked, nebParams, wave=wave_masked)

# Fit Stellar component
stellarMod = lmfit.models.ExponentialModel(prefix='stellar_')
stellarParams = stellarMod.guess(voxel_int_masked, x=wave_masked)
out_stellar = stellarMod.fit(voxel_int_masked, stellarParams, x=wave_masked)

# Fit Combined model
contMod = nebModel + stellarMod
nebParams.update(stellarParams)
# out_comb = contMod.fit(voxel_int_masked, nebParams, wave=wave_masked, x=wave_masked)
# comps = out_comb.eval_components(wave=wave_masked)
#
# # Plot the results
# fig, ax = plt.subplots(figsize=(12, 81))
# ax.plot(wave_masked, flux_masked * flux_norm, label='Continuum spectrum')
# ax.plot(wave_masked, out_stellar.best_fit, label='Stellar fit')
# ax.plot(wave_masked, out_neb.best_fit, label='Nebular fit')
# ax.plot(wave_masked, out_comb.best_fit, label='Combined fit')
# ax.plot(wave_masked, comps['neb_'], linestyle='--', label='Component neb')
# ax.plot(wave_masked, comps['stellar_'], linestyle='--', label='Component stellar')
# ax.update({'xlabel': r'Wavelength $(\AA)$', 'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})$'})
# ax.legend()
# plt.show()

# Plot the results
labelsDict = {'xlabel': r'Wavelength $(\AA)$',
              'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})\cdot10^{20}$',
              'title': f'Galaxy CGCG007 (voxel coords: {idx_voxel[0]}, {idx_voxel[1]})'}

fig, ax = plt.subplots(figsize=(12, 81))
ax.plot(wave, voxel_int, label='Spectrum intensity')
ax.plot(wave_masked, voxel_int_masked, label='Masked lines')
ax.plot(wave, neb_int, label='Nebular calculation')
# ax.plot(wave_masked, out_neb.best_fit, label='Nebular fit')
# ax.plot(wave_masked, out_stellar.best_fit, label='Stellar fit')

# ax.plot(wave_masked, out_comb.best_fit, label='Combined fit')
# ax.plot(wave_masked, comps['neb_'], linestyle='--', label='Component neb')
# ax.plot(wave_masked, comps['stellar_'], linestyle='--', label='Component stellar')
ax.update(labelsDict)
ax.set_yscale('log')
ax.legend()
plt.show()
