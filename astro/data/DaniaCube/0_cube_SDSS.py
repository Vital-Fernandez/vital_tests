import numpy as np
import lime

from pathlib import Path
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt

cube_path = Path(f'U:\AstroData\Dania_cube\J0021_4800_conOII_relativoSDSS.fits')
z_obj = 0.09835

fits.info(cube_path)

with fits.open(cube_path) as hdul:
    frame_idx = 0
    data, header = hdul[frame_idx].data, hdul[frame_idx].header

w_min = header['CRVAL3']
dw = header['CD3_3']  # dw = 0.862936 INDEF (Wavelength interval per pixel)
pixels = header['NAXIS3']  # nw = 3801 number of output pixels
w_max = w_min + dw * pixels
wave = np.linspace(w_min, w_max, pixels, endpoint=False)
flux = data
wcs = WCS(header)

J0021 = lime.Cube(wave, data, redshift=z_obj, wcs=wcs)
J0021.check.cube('H1_4861A', rest_frame=True)

spatial_mask_file = 'J0021_4800_spatial_mask.fits'
J0021.spatial_masker('O3_5007A', param='SN_line', percentiles=(80, 95),  output_address='J0021_4800_spatial_mask.fits')
J0021.plot.cube('O3_5007A', masks_file=spatial_mask_file)
J0021.check.cube('H1_4861A', rest_frame=True, masks_file=spatial_mask_file)

spaxel = J0021.get_spaxel(12, 6)
spaxel.fit.band('O3_5007A')
spaxel.plot.band()

log_measurements = 'J0021_4800_measurements.fits'
default_bands = lime.io._parent_bands_file
J0021.fit.spatial_mask(spatial_mask_file, bands_frame=default_bands, lines_list=['O3_5007A'], output_log=log_measurements)
J0021.check.cube('H1_4861A', lines_log_address=log_measurements, rest_frame=True)

# Export the log spaxel measurements as maps:
param_list = ['intg_flux', 'intg_err', 'gauss_flux', 'gauss_err']
lines_list = ['O3_5007A']
lime.save_parameter_maps(log_measurements, param_list, lines_list, output_folder='.', spatial_mask_file=spatial_mask_file,
                         output_files_prefix='J0021_', wcs=wcs)

# Recover integraetd flux map
intg_flux_map = fits.getdata('J0021_intg_flux.fits')
hdr = fits.getheader('J0021_intg_flux.fits', 'O3_5007A')
wcs_maps = WCS(hdr)

# Create the plot
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(projection=wcs_maps, slices=('x', 'y'))
im = ax.imshow(intg_flux_map)
cbar = fig.colorbar(im, ax=ax)
ax.update({'title': f'J0021 [OIII]5007A integrated flux', 'xlabel': r'RA', 'ylabel': r'DEC'})
plt.show()