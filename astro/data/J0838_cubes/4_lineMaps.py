import numpy as np
import src.specsiser as sr
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astropy.table import Table
import pyneb as pn
from src.specsiser.print.plot import STANDARD_PLOT

label_Conver = {'H1_6563A': 'Halpha',
                'H1_4861A': 'Hbeta',
                'H1_4341A': 'Hgamma',
                'H1_4102A': 'Hdelta'}

latex_Conver = {'H1_6563A': r'H\alpha',
                'H1_4861A': r'H\beta',
                'H1_4341A': r'H\gamma',
                'H1_4102A': r'H\delta',
                'v_r': r'$v_{r}\,(km/s)$',
                'sigma_vel': r'$\sigma_{int}\,(km/s)$'}

# Plot set up
defaultConf = STANDARD_PLOT.copy()
defaultConf['axes.titlesize'] = 20
rcParams.update(defaultConf)

# # Declare data and files location
# obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
# objList = obsData['data_location']['object_list']
# fileList = obsData['data_location']['file_list']
# fitsFolder = Path(obsData['data_location']['fits_folder'])
# dataFolder = Path(obsData['data_location']['data_folder'])
# resultsFolder = Path(obsData['data_location']['results_folder'])
# voxel_grid_size = obsData['sample_data']['grid_shape_array']
# coordinates_keys_list = obsData['data_location']['wcs_key_list']

# Declare data and files location
obsConf = sr.loadConfData('J0838_cubes.ini')
fitsFolder = Path(obsConf['data_location']['fits_folder'])
dataFolder = Path(obsConf['data_location']['data_folder'])
resultsFolder = Path(obsConf['data_location']['results_folder'])

fileList = obsConf['data_location']['file_list']
objList = obsConf['data_location']['object_list']
z_list = obsConf['sample_data']['z_array']
norm_flux = obsConf['sample_data']['norm_flux']
percentil_array = obsConf['sample_data']['percentil_array']
voxel_grid_size = obsConf['sample_data']['grid_shape_array']

# Store emissivity ratios at standard conditions
H1 = pn.RecAtom('H', 1)
temp, den = 10000.0, 100.0

theoEmis_dict = {}
for chemLabel, plotLabel in label_Conver.items():
    ion, wave, latexLabel = sr.label_decomposition(chemLabel, scalar_output=True)
    dict_label = f'{plotLabel}/Hdelta'
    theoRatio = H1.getEmissivity(temp, den, wave=wave) / H1.getEmissivity(temp, den, wave=4102)
    theoEmis_dict[dict_label] = theoRatio
red_model = sr.ExtinctionModel(Rv=obsConf['extinction']['R_v'], red_curve=obsConf['extinction']['red_law'])

# Data location
objFolder = resultsFolder
db_address = objFolder/'J0838_blue_database.fits'
linesMaps_fits_address = resultsFolder/f'J0838_lineParamMaps.fits'

# ----------------------------------------- Generate the image data
#
# # Empty containers for the images
# image_dict = {}
# for chemLabel, plotLabel in label_Conver.items():
#     image_dict[chemLabel] = np.full(voxel_grid_size.astype(int), np.nan)
# 
# # Fits output file
# new_hdul = fits.HDUList()
# new_hdul.append(fits.PrimaryHDU())
#
# hdr_plot = fits.getheader(db_address, extname='PlotConf')
# hdu_table = fits.BinTableHDU.from_columns(columns=[], header=hdr_plot, name='PlotConf')
# new_hdul.append(hdu_table)
# 
# for i, obj in enumerate(objList):
# 
#     # Data location
#     cube_address_i = fitsFolder / fileList[i]
#     mask_address = dataFolder / obsConf['data_location']['mask_global']
#     voxelFolder = resultsFolder / obj
#     color = 'blue' if 'blue' in obj else 'red'
#     fitsLog_addresss = objFolder / f'{obj}_lineslog.fits'
# 
#     # Extinction model
#     regions_to_treat = [0, 1, 2]
# 
#     # Open the lines log database
#     with fits.open(fitsLog_addresss) as hdul:
# 
#         # Loop throught the line regions
#         for idx_region in regions_to_treat:
#             region_label = f'region_{idx_region}'
#             region_mask = fits.getdata(db_address, region_label, ver=1)
#             region_mask = region_mask.astype(bool)
#             idcs_voxels = np.argwhere(region_mask)
# 
#             # Loop through the region voxels
#             for idx_voxel, idx_pair in enumerate(idcs_voxels):
#                 idx_j, idx_i = idx_pair
#                 logLabel = f'{idx_j}-{idx_i}_linelog'
# 
#                 # Load lines log data and store it as an image
#                 if logLabel in hdul:
#                     linesDF = Table(hdul[logLabel].data).to_pandas()
#                     linesDF.set_index('index', inplace=True)
# 
#                     for chemLabel, plotLabel in label_Conver.items():
#                         if chemLabel in linesDF.index:
#                             dict_label = f'{plotLabel}'
#                             lineFlux = linesDF.loc[chemLabel, 'intg_flux']
#                             # print(chemLabel, lineFlux)
#                             image_dict[chemLabel][idx_j, idx_i] = lineFlux# - theoEmis_dict[dict_label]
# 
# for param_line, param_map in image_dict.items():
#     new_hdul.append(fits.ImageHDU(name=param_line, data=param_map, ver=1))
# new_hdul.writeto(linesMaps_fits_address, overwrite=True, output_verify='fix')

# ----------------------------------------- Generate the image plots

hdr_plot = fits.getheader(linesMaps_fits_address, extname='PlotConf')
flux5007_image = fits.getdata(db_address, f'O3_5007A_FLUX', ver=1)
# linthresh = flux6563_levels[-2], vmin = flux6563_levels[-3],

rcParams.update(defaultConf)

halpha_cmap = cm.gray
hbeta_image = fits.getdata(linesMaps_fits_address, 'H1_4102A', ver=1)

# Recombination lines
for chemLabel, plotLabel in label_Conver.items():

    if chemLabel != 'H1_4102A':
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot()

        dict_label = f'{plotLabel}/Hdelta'
        flux_image = fits.getdata(linesMaps_fits_address, chemLabel, ver=1)/hbeta_image

        print(chemLabel)
        print(np.nanmin(flux_image), theoEmis_dict[dict_label], np.nanmax(flux_image))
        print()
        divnorm = colors.TwoSlopeNorm(vmin=np.nanmin(flux_image),
                                      vcenter=theoEmis_dict[dict_label],
                                      vmax=np.nanmax(flux_image))

        im = ax.imshow(flux5007_image, cmap=halpha_cmap)
        im2 = ax.imshow(flux_image, cmap='RdBu_r', norm=divnorm)

        cbar = fig.colorbar(im2, ax=ax)
        cbar.set_label('Line ratio (white theoretical value)', rotation=270, labelpad=50, fontsize=15)

        ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_Conver[chemLabel], latex_Conver['H1_4861A'])
        ax.update({'title': r'Galaxy J0838: {}'.format(ratio_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
        plt.tight_layout()
        plt.show()