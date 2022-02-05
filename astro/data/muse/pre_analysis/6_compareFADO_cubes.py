from pprint import pprint

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import src.specsiser as sr
from pathlib import Path
import astropy.units as u
from mpdaf.obj import deg2sexa
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
import src.specsiser as sr


def reconstruct_wavelength(header):
    dw = header['CDELT3']
    w_min = header['CRVAL3']
    nPixels = header['NAXIS3']
    w_max = w_min + dw * nPixels
    return np.linspace(w_min, w_max, nPixels, endpoint=False)



conf_file_address = '../muse_greenpeas.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)
fits_folder = Path('/home/vital/Astro-data/Observations/MUSE - Amorin/FADO_analysis/Z4SalpP2000/')
fits_1 = 'cgcg007025_HBIN024_FDres2_2DnebSED.fits'
fits_2 = 'cgcg007025_HBIN024_FDres2_2DnoNEB.fits'
fits_3 = 'cgcg007025_HBIN024_FDres2_3DnebSED.fits'
fits_4 = 'cgcg007025_HBIN024_FDres2_3DnoNEB.fits'
fits_5 = 'cgcg007025_HBIN024_FDres2_3DOBS.fits'
fits_6 = 'cgcg007025_HBIN024_FDres2_3DstelFIT.fits'

coord = (170, 170)

file_path = Path('/home/vital/Astro-data/Observations/MUSE - Amorin/CGCG007.fits')
wave, cube, header = sr.import_fits_data(file_path, instrument='MUSE')
z_obj = 0.004691

file_address, ext = fits_folder/fits_3, 0
with fits.open(file_address) as hdu_list:
    data3 = hdu_list[ext].data
    hdr3 = hdu_list[ext].header

file_address, ext = fits_folder/fits_4, 0
with fits.open(file_address) as hdu_list:
    data4 = hdu_list[ext].data
    hdr4 = hdu_list[ext].header

file_address, ext = fits_folder/fits_5, 0
with fits.open(file_address) as hdu_list:
    data5 = hdu_list[ext].data
    hdr5 = hdu_list[ext].header

file_address, ext = fits_folder/fits_6, 0
with fits.open(file_address) as hdu_list:
    data6 = hdu_list[ext].data
    hdr6 = hdu_list[ext].header

wave3 = reconstruct_wavelength(hdr3)
flux3 = data3[:, coord[0], coord[1]]

wave4 = reconstruct_wavelength(hdr4)
flux4 = data4[:, coord[0], coord[1]]

wave5 = reconstruct_wavelength(hdr5)
flux5 = data5[:, coord[0], coord[1]]

wave6 = reconstruct_wavelength(hdr6)
flux6 = data6[:, coord[0], coord[1]]


wave_voxel = wave / (1+z_obj)
flux_voxel = cube[:, coord[0], coord[1]].data.data * (1 + z_obj)


w_min, w_max = 4800, 8200
idcs_orig = np.searchsorted(wave_voxel, (w_min-1, w_max))
idcs_fado = np.searchsorted(wave5, (w_min, w_max))

wave_fado = wave5[idcs_fado[0]:idcs_fado[-1]]
flux_fado = flux5[idcs_fado[0]:idcs_fado[-1]]

orig_interpolator = interp1d(wave_voxel, flux_voxel)
flux_orig_interp = orig_interpolator(wave_fado)

# print(wave_voxel)
# print(wave5)

# flux_orig_full = orig_interpolator(wave5[5:-5])


const = flux_orig_interp/flux_fado
print(const.mean())
print(const.std())

fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(wave_voxel, flux_voxel, label=f'Original')
# ax.plot(wave3, flux3, label=f'3DnebSED')
# ax.plot(wave4, flux4, label=f'3DnoNEB')
ax.plot(wave5, flux5, label=f'3DOBS')
# ax.plot(wave6, flux6, label=f'3DstelFIT ')
ax.plot(wave_fado, flux_fado, label=f'3DOBS_section', linestyle='--')
ax.plot(wave_fado, flux_orig_interp/const.mean(), label=f'Orig interp', linestyle='--')

ax.set_yscale('log')
ax.legend()
ax.update({'xlabel': r'Wavelength $(\AA)$', r'ylabel': 'Fits normalized flux'})
plt.show()

# #------------- cgcg007025_HBIN024_FDres2_3DstelFIT
# file_address, ext = fits_folder/fits_6, 0
# with fits.open(file_address) as hdu_list:
#     data = hdu_list[ext].data
#     hdr = hdu_list[ext].header
# print(f'Treating file: {fits_6}')
# fits.info(file_address)
# print(f'- Shape: {data.shape}')
# print(f'- Header:')
# pprint(hdr)
#
# wave = reconstruct_wavelength(hdr)
# flux = data[:, coord[0], coord[1]]
#
# fig, ax = plt.subplots(figsize=(16, 8))
# ax.plot(wave, flux, label=f'Voxel flux ({fits_6})')
# ax.legend()
# ax.update({'xlabel': r'Wavelength (\AA)', 'ylabel': 'Flux', 'title': 'Gaussian fitting'})
# plt.show()