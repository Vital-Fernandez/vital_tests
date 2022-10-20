import numpy as np
from pathlib import Path
from astropy import units as u
import pysynphot as S
import matplotlib.pyplot as plt
import lime


def mJy_to_flam(wave, flux):

    wave_A = wave * 10
    flux_conv = flux * 2.99792458e-8 / np.square(wave_A)

    return flux_conv


wave_array = np.linspace(500, 900, 900-500)
flux_array = np.ones(wave_array.size)

s_spec = S.spectrum.ArraySourceSpectrum(wave_array, flux_array, waveunits='nm', fluxunits='mjy', keepneg=True)
print(f'{s_spec.waveunits.name} and {s_spec.fluxunits.name}')

# Converting the units
s_spec.convert('flam')
print(f'{s_spec.waveunits.name} and {s_spec.fluxunits.name}')
print(s_spec.flux[:4])

wave_astro = wave_array * u.nm
flux_astro = flux_array * u.mJy
flux_astro_flam = flux_astro.to(u.erg/u.s/u.cm**2/u.AA, equivalencies=u.spectral_density(wave_astro))
print(flux_astro_flam[:4])

flux_mio = mJy_to_flam(wave_array, flux_array)
print(flux_mio[:4])


# Converting to mJys
# def specPerAngs_to_specPerHz(spec_flux):
#     return
#
# # Original spectrum from xshooter
# fitsFileAddress = Path('D:/Downloads/ADP.2019-02-11T15_41_46.691.fits')
# wave_nm, flux_array, headers = lime.load_fits(fitsFileAddress, instrument='xshooter')
# wave_nm_pysyn, flux_array_pysyn, headers_pysyn = sr.import_fits_data(fitsFileAddress, instrument='xshooter')
#
# # Wavelength conversion to angstroms
# wave, flux = wave_nm * 10, flux_array[1]
# wave_pysyn, flux_pysyn =wave_nm_pysyn * 10, flux_array_pysyn[1]
#
# # Conversion from gemini page
# fitsFileAddress_gemini = Path('D:/Downloads/ADP.2019-02-11T15_41_46.691_conversion.txt')
# wave_gem, flux_gemn = np.loadtxt(fitsFileAddress_gemini, unpack=True)
#
# # ---- Using pysynphot
# spec = pys.spectrum.ArraySourceSpectrum(wave=wave_pysyn, flux=flux_pysyn, waveunits='angstroms', fluxunits='flam',
#                                         keepneg=True)
# spec.convert('mjy')
# spec.convert('micron')
#
# # Manual conversion
# wave_nm
# c_ang_per_s = 2.998e+18
# erg_2_mju = 1e26
# wave_um = 0.0001 * wave
# flux_mJy = (flux * np.power(wave, 2) / c_ang_per_s) * erg_2_mju
#
# # # Plot the input spectra
# # fig, ax = plt.subplots(figsize=(12, 8))
# # ax.step(wave, flux, label='Orig')
# # ax.step(wave_pysyn, flux_pysyn, label='Orig pysn')
# # ax.legend()
# # plt.tight_layout()
# # plt.show()
#
#
# # Plot the conversions
# fig, ax = plt.subplots(figsize=(12, 8))
# ax.step(spec.wave, spec.flux, label='spec pysn')
# # ax.step(wave_gem, flux_gemn, label='spec gemini')
# ax.step(wave_um, flux_mJy, label='my conversion', linestyle=':')
# ax.legend()
# plt.tight_layout()
# plt.show()
#
# # Save conversion
# output_address = 'D:/Downloads/ADP.2019-02-11T15_41_46.691_um-mJy.txt'
# np.savetxt(output_address, np.c_[spec.wave, spec.flux])
#
# # output_address = 'D:/Downloads/ADP.2019-02-11T15_41_46.691.txt'
# # np.savetxt(output_address, np.c_[wave_A2, flux_array2[1]])
#
#
# # unit = 1 * u.erg/(u.s * u.cm**2 * u.AA)
# # unit.to(u.erg/(u.s * u.cm**2 * u.Hz))
# #
# # c_ang_per_s = 2.998e+18
# # mjy = 10-26 erg cm-2 s-1 Hz-1
#
# # Pysynphot conversion
# # a = pys.spectrum.ArraySourceSpectrum(wave=wave_A, flux=flux_array[1], waveunits='angstroms', fluxunits='flam', keepneg=True)
# # a.convert('mjy')
# #
# # # Manual conversion
# # c_ang_per_s = 2.998e+18
# # flux_nu = flux_array[1] * wave_A2**2 / c_ang_per_s
# #
# #
# # lm = sr.LineMesurer(wave_A, flux_array[1], input_err=flux_array[2], redshift=0.301)
# # lm.plot_spectrum_components(continuumFlux=flux_array2[1])
