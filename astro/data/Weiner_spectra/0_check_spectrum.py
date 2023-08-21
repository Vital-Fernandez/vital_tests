import lime
import numpy as np

file_address = './data/pegase_e13_lum_norm.csv'

wave, flux1, flux2, flux3 = np.loadtxt(file_address, skiprows=1, unpack=True, delimiter=',')

print(wave)
for flux in [flux1, flux2, flux3]:
    spec = lime.Spectrum(wave, flux)
    spec.plot.spectrum(maximize=True)

file_address = './data/elg_from_sd_lum_norm.csv'

wave, flux1, flux2 = np.loadtxt(file_address, skiprows=1, unpack=True, delimiter=',')

print(wave)
spec = lime.Spectrum(wave, flux2)
spec.plot.spectrum(maximize=True)