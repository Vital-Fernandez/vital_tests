import lime
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt, cm, colors, patches
from astropy.wcs import WCS
from mpdaf.obj import Cube


# Data location
cube_address = f'/mnt/AstroData/Observations/Ana_cube/uboala.fits'
cfg_file = '/astro/collabs/Ana_Cube/config_musep.cfg'
mask_file = '/astro/collabs/Ana_Cube/mascara_musep.txt'
spatial_mask_file = f'/astro/collabs/Ana_Cube/mascara_musep.fits'
log_address = f'/home/vital/PycharmProjects/vital_tests/astro/data/Ana_Cube/line_measurements.fits'

# Load the configuration file:
obs_cfg = lime.load_cfg(cfg_file)
z_SHOC148 = 0.16481 # obs_cfg['SHOC148_data']['z']
norm_flux = 1e-20 # obs_cfg['SHOC148_data']['norm_flux']

# Percentiles for the spatial mask
percentile_array = np.array([85, 95, 99])

# Parent mask file for the lines:
mask_log = lime.load_lines_log(mask_file)

# Load MUSE cube
galaxy_IFU = Cube(cube_address)

# Crop the IFU to the galaxy location and remove the MUSE normalization
spec_zone = galaxy_IFU[:, 130:190, 130:190].data.data * norm_flux

# Calculate the wavelength array
header = galaxy_IFU.data_header
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

# Labels for the axes for the plot
ax_conf = {'image': {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'MUSE SHOC148'}}

# Color normalization for the flux band:
min_flux = np.nanpercentile(Halpha_image, 60)
log_norm_bg = colors.SymLogNorm(linthresh=min_flux, vmin=min_flux)

# Create a dictionary with the coordinate entries for the header
hdr_coords = {}
for key in lime.COORD_ENTRIES:
    if key in header:
        hdr_coords[key] = header[key]

# Plot of the Central voxel at 34 - 26, Interactive plotter for IFU data cubes
lime.CubeInspector(wave, spec_zone, Halpha_image, OIII_image, OIII_contourLevels, color_norm=log_norm_bg,
                   fits_header=header, ax_cfg=ax_conf)


# Generate the mask file
lime.spatial_mask_generator('SN_line', wave, spec_zone, percentile_array, signal_band=OIII_band,
                             cont_band=OIII_cont, mask_ref='O3_5007A', output_address=spatial_mask_file,
                             fits_header=hdr_coords, show_plot=True)

core_voxel_coord = (34, 26)
spec = lime.Spectrum(wave, spec_zone[:, core_voxel_coord[0], core_voxel_coord[1]], redshift=z_SHOC148,
                     norm_flux=norm_flux)
spec.plot_spectrum()
#
# # Loop throught the masks and analyse the data
# hdul_log = fits.HDUList([fits.PrimaryHDU()])
# for idx_mask in [0, 1, 2]:
#
#     mask_name = f'O3_5007A_MASK_{idx_mask}'
#
#     # Load the region spatial mask:
#     region_label = f'O3_5007A_MASK_{idx_mask}'
#     region_mask = fits.getdata(spatial_mask_file, region_label, ver=1)
#     region_mask = region_mask.astype(bool)
#
#     # Convert the mask into an array of spaxel coordinates (idxY, idxX)
#     idcs_spaxels = np.argwhere(region_mask)
#
#     # Load the region fitting configuration
#     region_fit_cfg = obs_cfg['line_fitting']
#
#     # Loop through the spaxels
#     print(f'- Treating region {idx_mask} with {np.sum(region_mask)} voxels')
#     for idx_spaxel, coords_spaxel in enumerate(idcs_spaxels):
#
#         # Define a spectrum object for the current spaxel
#         idxY, idxX = coords_spaxel
#         flux = spec_zone[:, idxY, idxX]
#         pixel_mask = np.isnan(flux)
#         spaxel_spec = lime.Spectrum(wave, spec_zone[:, idxY, idxX], redshift=z_SHOC148, norm_flux=norm_flux,
#                                     pixel_mask=pixel_mask)
#
#         # Limit the line fittings to those detected
#         matched_mask_log = spaxel_spec.line_detection(lines_log=mask_log)
#
#         # Loop through the detected lines
#         print(f'-- Treating spaxel {idx_spaxel}')
#         for idx_line, line in enumerate(matched_mask_log.index):
#
#             wave_regions = matched_mask_log.loc[line, 'w1':'w6'].values
#             try:
#                 spaxel_spec.fit_from_wavelengths(line, wave_regions, fit_method='least_squares', user_cfg=region_fit_cfg)
#
#             except ValueError as e:
#                 print(f'--- Line measuring failure at {line} in spaxel {idxY}-{idxX}:\n{e}')
#
#         spaxel_spec.plot_spectrum(include_fits=True)
#
#         # Convert the measurements log into a HDU and append it to the HDU list unless it is empty
#         linesHDU = lime.log_to_HDU(spaxel_spec.log, ext_name=f'{idxY}-{idxX}_LINESLOG', header_dict=hdr_coords)
#
#         # Check the HDU is not empty (no lines measured)
#         if linesHDU is not None:
#             hdul_log.append(linesHDU)
#
#     # After the regions spaxels have been analysed save all the measurements to a .fits file
#     hdul_log.writeto(log_address, overwrite=True, output_verify='fix')
