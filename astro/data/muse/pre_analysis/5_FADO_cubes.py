from pprint import pprint

import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
import astropy.units as u
from mpdaf.obj import deg2sexa
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt


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

#------------- HBIN024_FDres2_2DnebSED
file_address, ext = fits_folder/fits_1, 0
with fits.open(file_address) as hdu_list:
    data = hdu_list[ext].data
    hdr = hdu_list[ext].header
print(f'Treating file: {fits_1}')
fits.info(file_address)
print(f'- Shape: {data.shape}')
print(f'- Header:')
pprint(hdr)
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y'))
im = ax.imshow(data)
cbar = fig.colorbar(im)
ax.update({'title': f'{fits_1}'})
plt.show()

#------------- cgcg007025_HBIN024_FDres2_2DnoNEB
file_address, ext = fits_folder/fits_2, 0
with fits.open(file_address) as hdu_list:
    data = hdu_list[ext].data
    hdr = hdu_list[ext].header
print(f'Treating file: {fits_2}')
fits.info(file_address)
print(f'- Shape: {data.shape}')
print(f'- Header:')
pprint(hdr)
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y'))
im = ax.imshow(data)
cbar = fig.colorbar(im)
ax.update({'title': f'{fits_2}'})
plt.show()

#------------- cgcg007025_HBIN024_FDres2_3DnebSED
file_address, ext = fits_folder/fits_3, 0
with fits.open(file_address) as hdu_list:
    data = hdu_list[ext].data
    hdr = hdu_list[ext].header
print(f'Treating file: {fits_3}')
fits.info(file_address)
print(f'- Shape: {data.shape}')
print(f'- Header:')
pprint(hdr)

wave = reconstruct_wavelength(hdr)
flux = data[:, coord[0], coord[1]]

fig, ax = plt.subplots(figsize=(16,8))
ax.plot(wave, flux, label=f'Voxel flux ({fits_3})')
ax.legend()
ax.update({'xlabel': r'Wavelength (\AA)', 'ylabel': 'Flux', 'title': 'Gaussian fitting'})
plt.show()

#------------- cgcg007025_HBIN024_FDres2_3DnoNEB
file_address, ext = fits_folder/fits_4, 0
with fits.open(file_address) as hdu_list:
    data = hdu_list[ext].data
    hdr = hdu_list[ext].header
print(f'Treating file: {fits_4}')
fits.info(file_address)
print(f'- Shape: {data.shape}')
print(f'- Header:')
pprint(hdr)

wave = reconstruct_wavelength(hdr)
flux = data[:, coord[0], coord[1]]

fig, ax = plt.subplots(figsize=(16,8))
ax.plot(wave, flux, label=f'Voxel flux ({fits_4})')
ax.legend()
ax.update({'xlabel': r'Wavelength (\AA)', 'ylabel': 'Flux', 'title': 'Gaussian fitting'})
plt.show()

#------------- cgcg007025_HBIN024_FDres2_3DOBS
file_address, ext = fits_folder/fits_5, 0
with fits.open(file_address) as hdu_list:
    data = hdu_list[ext].data
    hdr = hdu_list[ext].header
print(f'Treating file: {fits_5}')
fits.info(file_address)
print(f'- Shape: {data.shape}')
print(f'- Header:')
pprint(hdr)

wave = reconstruct_wavelength(hdr)
flux = data[:, coord[0], coord[1]]

fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(wave, flux, label=f'Voxel flux ({fits_5})')
ax.legend()
ax.update({'xlabel': r'Wavelength (\AA)', 'ylabel': 'Flux', 'title': 'Gaussian fitting'})
plt.show()

#------------- cgcg007025_HBIN024_FDres2_3DstelFIT
file_address, ext = fits_folder/fits_6, 0
with fits.open(file_address) as hdu_list:
    data = hdu_list[ext].data
    hdr = hdu_list[ext].header
print(f'Treating file: {fits_6}')
fits.info(file_address)
print(f'- Shape: {data.shape}')
print(f'- Header:')
pprint(hdr)

wave = reconstruct_wavelength(hdr)
flux = data[:, coord[0], coord[1]]

fig, ax = plt.subplots(figsize=(16, 8))
ax.plot(wave, flux, label=f'Voxel flux ({fits_6})')
ax.legend()
ax.update({'xlabel': r'Wavelength (\AA)', 'ylabel': 'Flux', 'title': 'Gaussian fitting'})
plt.show()