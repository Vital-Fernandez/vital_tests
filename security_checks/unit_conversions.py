import numpy as np
import src.specsiser as sr
from pathlib import Path
from astropy import units as u
import pysynphot as pys
import matplotlib.pyplot as plt

def specPerAngs_to_specPerHz(spec_flux):
    return

# Original spectrum from xshooter
fitsFileAddress = Path('D:/Downloads/ADP.2019-02-11T15_41_46.691.fits')
wave_nm, flux_array, headers = sr.import_fits_data(fitsFileAddress, instrument='xshooter')
wave_nm_pysyn, flux_array_pysyn, headers_pysyn = sr.import_fits_data(fitsFileAddress, instrument='xshooter')

# Wavelength conversion to angstroms
wave, flux = wave_nm * 10, flux_array[1]
wave_pysyn, flux_pysyn =wave_nm_pysyn * 10, flux_array_pysyn[1]

# Conversion from gemini page
fitsFileAddress_gemini = Path('D:/Downloads/ADP.2019-02-11T15_41_46.691_conversion.txt')
wave_gem, flux_gemn = np.loadtxt(fitsFileAddress_gemini, unpack=True)

# ---- Using pysynphot
spec = pys.spectrum.ArraySourceSpectrum(wave=wave_pysyn, flux=flux_pysyn, waveunits='angstroms', fluxunits='flam',
                                        keepneg=True)
spec.convert('mjy')
spec.convert('micron')

# Manual conversion
wave_nm
c_ang_per_s = 2.998e+18
erg_2_mju = 1e26
wave_um = 0.0001 * wave
flux_mJy = (flux * np.power(wave, 2) / c_ang_per_s) * erg_2_mju

# # Plot the input spectra
# fig, ax = plt.subplots(figsize=(12, 8))
# ax.step(wave, flux, label='Orig')
# ax.step(wave_pysyn, flux_pysyn, label='Orig pysn')
# ax.legend()
# plt.tight_layout()
# plt.show()


# Plot the conversions
fig, ax = plt.subplots(figsize=(12, 8))
ax.step(spec.wave, spec.flux, label='spec pysn')
# ax.step(wave_gem, flux_gemn, label='spec gemini')
ax.step(wave_um, flux_mJy, label='my conversion', linestyle=':')
ax.legend()
plt.tight_layout()
plt.show()

# Save conversion
output_address = 'D:/Downloads/ADP.2019-02-11T15_41_46.691_um-mJy.txt'
np.savetxt(output_address, np.c_[spec.wave, spec.flux])

# output_address = 'D:/Downloads/ADP.2019-02-11T15_41_46.691.txt'
# np.savetxt(output_address, np.c_[wave_A2, flux_array2[1]])


# unit = 1 * u.erg/(u.s * u.cm**2 * u.AA)
# unit.to(u.erg/(u.s * u.cm**2 * u.Hz))
#
# c_ang_per_s = 2.998e+18
# mjy = 10-26 erg cm-2 s-1 Hz-1

# Pysynphot conversion
# a = pys.spectrum.ArraySourceSpectrum(wave=wave_A, flux=flux_array[1], waveunits='angstroms', fluxunits='flam', keepneg=True)
# a.convert('mjy')
#
# # Manual conversion
# c_ang_per_s = 2.998e+18
# flux_nu = flux_array[1] * wave_A2**2 / c_ang_per_s
#
#
# lm = sr.LineMesurer(wave_A, flux_array[1], input_err=flux_array[2], redshift=0.301)
# lm.plot_spectrum_components(continuumFlux=flux_array2[1])
