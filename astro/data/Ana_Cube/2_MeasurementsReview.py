import lime
import numpy as np
from astropy.io import fits
from matplotlib import colors
from astropy.wcs import WCS
from mpdaf.obj import Cube


# Data location
cube_address = f'/media/vital/Vital stick/Ana/SHOC148_MUSE_cube.fits'
cfg_file = '/home/vital/PycharmProjects/vital_tests/astro/data/Ana_Cube/config_musep.cfg'
mask_file = '/home/vital/PycharmProjects/vital_tests/astro/data/Ana_Cube/mascara_musep.txt'
spatial_mask_file = f'/home/vital/PycharmProjects/vital_tests/astro/data/Ana_Cube/mascara_musep.fits'
log_address = f'/home/vital/PycharmProjects/vital_tests/astro/data/Ana_Cube/line_measurements.fits'

# Load the configuration file:
obs_cfg = lime.load_cfg(cfg_file)
z_SHOC148 = obs_cfg['SHOC148_data']['z']
norm_flux = obs_cfg['SHOC148_data']['norm_flux']

# Percentiles para la mascara spacial
percentile_array = np.array([85, 96, 99])

# and the masks file:
mask_log = lime.load_lines_log(mask_file)

# Load MUSE cube
alfa = Cube(cube_address)

# Zona extendida
spec_zone = alfa[:, 130:190, 130:190].data.data * norm_flux

# Computing wavelength array
header = alfa.data_header
dw = header['CD3_3']
w_min = header['CRVAL3']
nPixels = header['NAXIS3']
w_max = w_min + dw * nPixels
wave = np.linspace(w_min, w_max, nPixels, endpoint=False)

# Using the Halpha band for the background image
Halpha_band = mask_log.loc['H1_6563A_b', 'w3':'w4'].values * (1 + z_SHOC148)
idcs_Halpha = np.searchsorted(wave, Halpha_band)
Halpha_image = spec_zone[idcs_Halpha[0]:idcs_Halpha[1], :, :].sum(axis=0)
Halpha_contourLevels = np.nanpercentile(Halpha_image, percentile_array)

# Use OIII lines as the foreground image contours
OIII_band = mask_log.loc['O3_5007A', 'w3':'w4'].values * (1 + z_SHOC148)
OIII_cont = mask_log.loc['O3_5007A', 'w1':'w2'].values * (1 + z_SHOC148)
idcs_OIII = np.searchsorted(wave, OIII_band)
OIII_image = spec_zone[idcs_OIII[0]:idcs_OIII[1], :, :].sum(axis=0)
OIII_contourLevels = np.nanpercentile(OIII_image, percentile_array)

# Labels for the axes
ax_conf = {'image': {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'MUSE SHOC148'}}

# Color normalization for the flux band:
min_flux = np.nanpercentile(Halpha_image, 60)
log_norm_bg = colors.SymLogNorm(linthresh=min_flux, vmin=min_flux)

# Create a dictionary with the coordinate entries for the header
hdr_coords = {}
for key in lime.COORD_ENTRIES:
    if key in header:
        hdr_coords[key] = header[key]

# Central voxel 34 - 26
# Interactive plotter for IFU data cubes
lime.CubeInspector(wave, spec_zone, Halpha_image, OIII_image, OIII_contourLevels, color_norm=log_norm_bg,
                   fits_header=header, ax_cfg=ax_conf, lines_log_address=log_address, redshift=z_SHOC148 )

