import numpy as np
import lime
from astropy.io import fits
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS

file_address = Path(r'S:\Astro_data\Observations\Dania_cube\J0926_R831_29julio.fits')

frame_idx = 0
with fits.open(file_address) as hdul:
    data, header = hdul[frame_idx].data, hdul[frame_idx].header

w_min = header['CRVAL3']
dw = header['CD3_3']  # dw = 0.862936 INDEF (Wavelength interval per pixel)
pixels = header['NAXIS3']  # nw = 3801 number of output pixels
w_max = w_min + dw * pixels
wave = np.linspace(w_min, w_max, pixels, endpoint=False)

print(data.shape)
region_Halpha = np.array([7732.0, 7764.0])

idcsHalpha = np.searchsorted(wave, region_Halpha)

wave_Halpha, fluxHalpha = wave[idcsHalpha[0]:idcsHalpha[1]], data[idcsHalpha[0]:idcsHalpha[1], :, :]

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot()
im = ax.imshow(fluxHalpha.sum(axis=1), aspect=0.5)
print(fluxHalpha.sum(axis=1).shape)

ax.update({'title': r'J0926 galaxy', 'xlabel': r'X', 'ylabel': r'Y'})
plt.show()

j_idx, i_idx = 10, 10
voxel_spec = lime.Spectrum(wave, data[:, j_idx, i_idx])
voxel_spec.plot_spectrum()

ext = f'coord-{coord_y}-{coord_x}'