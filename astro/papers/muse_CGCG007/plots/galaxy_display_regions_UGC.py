import lime
import numpy as np
import pandas as pd
import lime as lm
from pathlib import Path
from astro.data.muse.common_methods import lineAreas, store_frame_to_fits
from astropy.io import fits
from matplotlib import pyplot as plt, cm, colors, patches, rcParams
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits
from lime.plots import STANDARD_PLOT


def image_from_mask(mask_list, bg_value=0):

    image_array = np.full(mask_list[0].shape, bg_value)

    image_mask_levels = np.zeros(len(mask_list))
    for i, mask_pixels_array in enumerate(mask_list):
        print(np.sum(mask_pixels_array))
        image_array[mask_pixels_array] = bg_value + i
        image_mask_levels[i] = bg_value + i

    return image_array, image_mask_levels


# Declare data and files location
obsData = lm.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
plotsFolder = Path(obsData['data_location']['plots_folder'])
z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

dict_errs = {}
dict_nan_values = {}

ref_flux_line = 'S3_6312A'

arrow_style = dict(arrowstyle='->', color='yellow', connectionstyle='arc3,rad=-0.5', alpha=0.5)

# for i, obj in enumerate(objList):
#
#     # Data location
#     cube_address = fitsFolder/fileList[i]
#     objFolder = resultsFolder/obj
#     voxelFolder = resultsFolder/obj/'voxel_data'
#     db_addresss = objFolder/f'{obj}_database.fits'
#     mask_address = dataFolder/f'{obj}_mask.txt'
#
#     # Output data:
#     masks_plot = plotsFolder/f'{obj}_regions.png'
#
#     # Load data
#     wave, cube, header = import_muse_fits(cube_address)
#     wave_rest = wave / (1 + z_objs[i])
#     mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)
#
#     # Declare voxels to analyse
#     invers_pertil_array = pertil_array[::-1]
#     flux6312_image = fits.getdata(db_addresss, f'{ref_flux_line}_flux', ver=1)
#     flux6312_levels = np.nanpercentile(flux6312_image, invers_pertil_array)
#
#     flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
#     flux6563_levels = np.nanpercentile(flux6563_image, invers_pertil_array)
#
#     Halpha_idxMinPercentil = 6
#     Halpha_min_level = flux6563_levels[5]
#     SIII_idxMinPercentil = 3
#
#     # Plot combined mask
#     defaultConf = STANDARD_PLOT.copy()
#     rcParams.update(defaultConf)
#
#     # fig = plt.figure(figsize=(12, 8), dpi=600)
#     fig = plt.figure(figsize=(12, 8))
#     ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
#     im = ax.imshow(flux6563_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=flux6563_levels[-2],
#                                                                         vmin=flux6563_levels[-2],
#                                                                         base=10))
#     ax.update({'xlabel': r'RA', 'ylabel': r'DEC'})
#     ax.set_xlim(90, 220)
#     ax.set_ylim(70, 240)
#     plt.show()
#     # plt.savefig(masks_plot)




# Data location
cube_address = '/mnt/AstroData/Observations/Muse - Amorin/UGC5205.fits'

# Load data
z_obj = 0.005113
wave, cube, header = import_muse_fits(cube_address)
wave_rest = wave / (1 + z_obj)
mask_df = lime.spectral_mask_generator()

# # Declare voxels to analyse
# invers_pertil_array = pertil_array[::-1]
# flux6312_image = fits.getdata(db_addresss, f'{ref_flux_line}_flux', ver=1)
# flux6312_levels = np.nanpercentile(flux6312_image, invers_pertil_array)
#
# flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
# flux6563_levels = np.nanpercentile(flux6563_image, invers_pertil_array)

lineLimits = mask_df.loc['H1_6563A', 'w3':'w4'].values
line_image = cube.get_image(lineLimits * (1 +z_obj), subtract_off=True)
flux_image = line_image.data.data
invers_pertil_array = pertil_array[::-1]
levelContours = np.nanpercentile(flux_image, invers_pertil_array)

Halpha_min_level = levelContours[5]

# Plot combined mask
defaultConf = STANDARD_PLOT.copy()
rcParams.update(defaultConf)

# fig = plt.figure(figsize=(12, 8), dpi=600)
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
im = ax.imshow(flux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=levelContours[-21],
                                                                vmin=levelContours[-1],
                                                                base=10))
ax.update({'xlabel': r'RA', 'ylabel': r'DEC'})
# ax.set_xlim(90, 220)
# ax.set_ylim(70, 240)
plt.show()