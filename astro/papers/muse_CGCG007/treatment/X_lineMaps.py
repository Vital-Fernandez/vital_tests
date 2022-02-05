import lime
import numpy as np
import pyneb as pn
import lime as lm

from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astropy.table import Table

from src.specsiser.print.plot import STANDARD_PLOT
from astro.papers.muse_CGCG007.muse_CGCG007_methods import label_Conver, dinamicLines, import_muse_fits, latex_Conver
from src.specsiser.physical_model.extinction_model import ExtinctionModel


import sys

# def printProgressBar(i, max, postText, n_bar=10):
#
#     #size of progress bar
#     j = i/max
#     sys.stdout.write('\r')
#     message = f"[{'=' * int(n_bar * j):{n_bar}s}] {int(100 * j)}%  {postText}"
#     sys.stdout.write(message)
#     sys.stdout.flush()
#
#     return
#
# def parameter_maps(log_file_address, param_dict, output_folder, mask_file_address=None, ext_mask='all', image_shape=None,
#                    ext_log='_LINESLOG', default_value=np.nan, fits_prefix='', page_hdr={}):
#
#     assert Path(log_file_address).is_file(), f'- ERROR: lines log at {log_file_address} not found'
#     assert Path(output_folder).is_dir(), f'- ERROR: Output parameter maps folder {output_folder} not found'
#
#     # Compile the list of voxels to recover the provided masks
#     if mask_file_address is not None:
#
#         assert Path(mask_file_address).is_file(), f'- ERROR: mask file at {mask_file_address} not found'
#
#         with fits.open(mask_file_address) as maskHDUs:
#
#             # Get the list of mask extensions
#             if ext_mask == 'all':
#                 if ('PRIMARY' in maskHDUs) and (len(maskHDUs) > 1):
#                     mask_list = []
#                     for i, HDU in enumerate(maskHDUs):
#                         mask_name = HDU.name
#                         if mask_name != 'PRIMARY':
#                             mask_list.append(mask_name)
#                     mask_list = np.array(mask_list)
#                 else:
#                     mask_list = np.array(['PRIMARY'])
#             else:
#                 mask_list = np.array(ext_mask, ndmin=1)
#
#             # Combine all the mask voxels into one
#             for i, mask_name in enumerate(mask_list):
#                 if i == 0:
#                     mask_array = maskHDUs[mask_name].data
#                     image_shape = mask_array.shape
#                 else:
#                     assert image_shape == maskHDUs[mask_name].data.shape, '- ERROR: Input masks do not have the same dimensions'
#                     mask_array += maskHDUs[mask_name].data
#
#             # Convert to boolean
#             mask_array = mask_array.astype(bool)
#
#             # List of spaxels in list [(idx_j, idx_i), ...] format
#             spaxel_list = np.argwhere(mask_array)
#
#     # No mask file is provided and the user just defines an image size tupple (nY, nX)
#     else:
#         mask_array = np.ones(image_shape).astype(bool)
#         spaxel_list = np.argwhere(mask_array)
#
#     # Generate containers for the data:
#     images_dict, param_list = {}, []
#     for param, line_list in param_dict.items():
#
#         # Make sure is an array and loop throuh them
#         line_list = np.array(line_list, ndmin=1)
#         for line in line_list:
#             images_dict[f'{param}-{line}'] = np.full(image_shape, default_value)
#             param_list.append([param, line])
#
#     param_list = np.array(param_list)
#
#     # Loop through the spaxels and fill the parameter images
#     n_spaxels = spaxel_list.shape[0]
#     spaxel_range = np.arange(n_spaxels)
#
#     bench_array = np.array([0, 25, 50, 75])
#     percentil_bench = np.percentile(spaxel_range, bench_array).astype(int)
#     with fits.open(log_file_address) as logHDUs:
#
#         for i_spaxel in spaxel_range:
#             idx_j, idx_i = spaxel_list[i_spaxel]
#             spaxel_ref = f'{idx_j}-{idx_i}{ext_log}'
#
#             printProgressBar(i_spaxel, n_spaxels, postText=f'spaxels treated ({n_spaxels})')
#
#             # Confirm log extension exists
#             if spaxel_ref in logHDUs:
#
#                 # Recover extension data
#                 log_data = logHDUs[spaxel_ref].data
#                 log_lines = log_data['index']
#
#                 # Loop through the parameters and the lines:
#                 for param, user_lines in param_dict.items():
#                     idcs_log = np.argwhere(np.in1d(log_lines, user_lines))
#                     for i_line in idcs_log:
#                         images_dict[f'{param}-{log_lines[i_line][0]}'][idx_j, idx_i] = log_data[param][i_line][0]
#
#     # New line after the rustic progress bar
#     print()
#
#     # Save the parameter maps as individual fits files with one line per page
#     for param, user_lines in param_dict.items():
#
#         # Primary header
#         paramHDUs = fits.HDUList()
#         paramHDUs.append(fits.PrimaryHDU())
#
#         # ImageHDU for the parameter maps
#         for line in user_lines:
#             hdr = fits.Header({'PARAM': param, 'LINE': param})
#             hdr.update(page_hdr)
#             data = images_dict[f'{param}-{line}']
#             paramHDUs.append(fits.ImageHDU(name=line, data=data, header=hdr, ver=1))
#
#         # Write to new file
#         output_file = Path(output_folder)/f'{fits_prefix}{param}.fits'
#         paramHDUs.writeto(output_file, overwrite=True, output_verify='fix')
#
#     return

# Plot set up
defaultConf = STANDARD_PLOT.copy()
defaultConf['axes.titlesize'] = 20
rcParams.update(defaultConf)

# Declare data and files location
obsData = lm.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
voxel_grid_size = obsData['sample_data']['grid_shape_array']
coordinates_keys_list = obsData['data_location']['wcs_key_list']

# Store emissivity ratios at standard conditions
H1 = pn.RecAtom('H', 1)
temp, den = 10000.0, 100.0
background_color = np.array((43, 43, 43))/255.0

theoEmis_dict = {}
for chemLabel, plotLabel in label_Conver.items():
    ion, wave, latexLabel = lm.label_decomposition(chemLabel, scalar_output=True)
    dict_label = f'{plotLabel}/Hbeta'
    theoRatio = H1.getEmissivity(temp, den, wave=wave) / H1.getEmissivity(temp, den, wave=4861)
    theoEmis_dict[dict_label] = theoRatio

# regions_to_treat = [0, 1]
regions_to_treat = [0, 1, 2, 3, 4, 5]

for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_address = objFolder / f'{obj}_database.fits'
    fitsLog_address = objFolder / f'{obj}_linesLog.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Load data
    wave, cube, header = import_muse_fits(cube_address)

    # Plot data
    hdr_plot = fits.getheader(db_address, extname='PlotConf')

    plot_dict = {'CRPIX1': header['CRPIX1'],
                 'CRPIX2': header['CRPIX2'],
                 'CD1_1': header['CD1_1'],
                 'CD1_2': header['CD1_2'],
                 'CD2_1': header['CD2_1'],
                 'CD2_2': header['CD2_2'],
                 'CUNIT1': header['CUNIT1'],
                 'CUNIT2': header['CUNIT2'],
                 'CTYPE1': header['CTYPE1'],
                 'CTYPE2': header['CTYPE2']}

    # Output file
    linesMaps_fits_address = objFolder/f'{obj}_lineParamMaps.fits'

    # Extinction model
    red_model = ExtinctionModel(Rv=obsData['Extinction']['R_v'], red_curve=obsData['Extinction']['red_law'])

    # Generate the image data
    param_images = {'intg_flux': ['H1_4861A', 'H1_6563A'],
                    'gauss_flux': ['H1_4861A', 'H1_6563A']}

    masks = ['MASK_0', 'MASK_1', 'MASK_2']
    lime.save_param_maps(fitsLog_address, param_images, objFolder, maskFits_address, ext_mask=masks, ext_log='_LINELOG',
                         page_hdr=plot_dict)

    # parameter_maps(fitsLog_address, param_images, objFolder, ext_log='_LINELOG',
    #                page_hdr=plot_dict, image_shape=(500,500))

    for param, user_lines in param_images.items():
        fits_file = Path(objFolder) / f'{param}.fits'
        with fits.open(fits_file):
            for line in user_lines:
                param_image = fits.getdata(fits_file, line)
                param_hdr = fits.getheader(fits_file, line)
                fig = plt.figure(figsize=(5, 5))
                ax = fig.add_subplot(projection=WCS(fits.Header(param_hdr)), slices=('x', 'y'))
                im = ax.imshow(param_image)
                ax.update({'title': f'Galaxy {obj}: {param}-{line}', 'xlabel': r'RA', 'ylabel': r'DEC'})
                # plt.tight_layout()
                plt.show()



    # # ----------------------------------------- Generate the image data
    #
    # # Empty containers for the images
    # image_dict = {}
    # for chemLabel, plotLabel in label_Conver.items():
    #     dict_label = f'{plotLabel}/Hbeta'
    #     image_dict[dict_label] = np.full(voxel_grid_size.astype(int), np.nan)
    #
    # for dinLabel, plotLabel in dinamicLines.items():
    #     for param in ('v_r', 'sigma_vel'):
    #         image_dict[f'{param}_{dinLabel}'] = np.full(voxel_grid_size.astype(int), np.nan)
    #
    # # Open the lines log database
    # with fits.open(fitsLog_address) as hdul:
    #
    #     # Loop throught the line regions
    #     for idx_region in regions_to_treat:
    #         region_label = f'region_{idx_region}'
    #         region_mask = fits.getdata(maskFits_address, region_label, ver=1)
    #         region_mask = region_mask.astype(bool)
    #         idcs_voxels = np.argwhere(region_mask)
    #
    #         # Loop through the region voxels
    #         for idx_voxel, idx_pair in enumerate(idcs_voxels):
    #             idx_j, idx_i = idx_pair
    #             logLabel = f'{idx_j}-{idx_i}_linelog'
    #
    #             # Load lines log data and store it as an image
    #             if logLabel in hdul:
    #                 linesDF = Table(hdul[logLabel].data).to_pandas()
    #                 linesDF.set_index('index', inplace=True)
    #
    #                 if 'H1_4861A' in linesDF.index:
    #                     Hbeta_flux = linesDF.loc['H1_4861A', 'gauss_flux']
    #                     for chemLabel, plotLabel in label_Conver.items():
    #
    #                         if chemLabel in linesDF.index:
    #                             dict_label = f'{plotLabel}/Hbeta'
    #                             lineFlux = linesDF.loc[chemLabel, 'gauss_flux']/Hbeta_flux
    #                             image_dict[dict_label][idx_j, idx_i] = lineFlux# - theoEmis_dict[dict_label]
    #
    #                 for dinLabel, plotLabel in dinamicLines.items():
    #                     if dinLabel in linesDF.index:
    #                         for param in ('v_r', 'sigma_vel'):
    #                             image_dict[f'{param}_{dinLabel}'][idx_j, idx_i] = linesDF.loc[dinLabel, param]
    #
    # # Storing the dictionary as a fits image file
    # new_hdul = fits.HDUList()
    # new_hdul.append(fits.PrimaryHDU())
    #
    # # Second page for the fits file plot configuration
    # col_waves = fits.Column(name='wave', array=wave, format='1E')
    # hdu_table = fits.BinTableHDU.from_columns([col_waves], name='PlotConf')
    # new_hdul.append(hdu_table)
    # for key in coordinates_keys_list:
    #     new_hdul[1].header[key] = cube.data_header[key]
    # new_hdul[1].header['NPIXWAVE'] = cube.data_header['NAXIS3']
    #
    # # Remaining pages for the line maps
    # for param_line, param_map in image_dict.items():
    #     new_hdul.append(fits.ImageHDU(name=param_line, data=param_map, ver=1))
    # new_hdul.writeto(linesMaps_fits_address, overwrite=True, output_verify='fix')


    # ----------------------------------------- Generate the image plots
    # hdr_plot = fits.getheader(linesMaps_fits_address, extname='PlotConf')
    # flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    # halpha_min_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')
    # halpha_thresd_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    # # linthresh = flux6563_levels[-2], vmin = flux6563_levels[-3],
    #
    # # defaultConf = DARK_PLOT.copy()
    # # rcParams.update(defaultConf)
    #
    # halpha_cmap = cm.gray
    # halpha_cmap.set_under(background_color)
    #
    # # Recombination lines
    # for chemLabel, plotLabel in label_Conver.items():
    #
    #     fig = plt.figure(figsize=(10, 10))
    #     ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))
    #
    #     dict_label = f'{plotLabel}/Hbeta'
    #     flux_image = fits.getdata(linesMaps_fits_address, dict_label, ver=1)
    #
    #     # divnorm = colors.TwoSlopeNorm(vmin=np.nanmin(flux_image),
    #     #                               vcenter=theoEmis_dict[dict_label],
    #     #                               vmax=np.nanmax(flux_image))
    #
    #     im = ax.imshow(flux6563_image, cmap=halpha_cmap,
    #                    norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))
    #     im2 = ax.imshow(flux_image, cmap='RdBu')#, norm=divnorm)
    #
    #     cbar = fig.colorbar(im2, ax=ax)
    #     cbar.set_label('Line ratio (white theoretical value)', rotation=270, labelpad=50, fontsize=15)
    #
    #     # ion, wavelength, latexLabel = lm.label_decomposition(latex_Conver[chemLabel], scalar_output=True)
    #
    #     ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_Conver[chemLabel], latex_Conver['H1_4861A'])
    #     ax.update({'title': r'Galaxy {}: {}'.format(obj, ratio_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    #     plt.tight_layout()
    #     # plt.show()
    #     plt.savefig(resultsFolder/obj/f'map_{obj}_{chemLabel}.png')
    #
    # # Line kinematics
    # vr_Halpha = fits.getdata(linesMaps_fits_address, 'v_r_H1_6563A', ver=1)
    # Halpha_label = dinamicLines['H1_6563A']
    # Halpha_mean, Halpha_std = np.nanmean(vr_Halpha), np.nanstd(vr_Halpha)
    #
    # print(Halpha_mean, Halpha_std)
    # for dinLabel, latex_label in dinamicLines.items():
    #     for param in ('v_r', 'sigma_vel'):
    #
    #         fig = plt.figure(figsize=(10, 10))
    #
    #         ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))
    #
    #         dict_label = f'{param}_{dinLabel}'
    #         param_image = fits.getdata(linesMaps_fits_address, dict_label, ver=1)
    #
    #         im = ax.imshow(flux6563_image, cmap=halpha_cmap,
    #                        norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))
    #
    #         if param == 'v_r':
    #             param_image = param_image - Halpha_mean
    #             param_min, param_max = np.nanmin(param_image), np.nanmax(param_image)
    #             divnorm = colors.TwoSlopeNorm(vcenter=0.0, vmin=-30, vmax=30)
    #             im2 = ax.imshow(param_image, cmap='RdBu_r', norm=divnorm)
    #             # title_label = f'$v_{{r}}$ {latex_label} - {Halpha_label}'
    #             title_label = f'$v_{{r}}$ {latex_label} -' + r'$\overline{v_{H\alpha}}$' + r' $({:.0f} \pm {:.0f}\,km/s)$'.format(Halpha_mean, Halpha_std)
    #
    #         else:
    #             im2 = ax.imshow(param_image)
    #             title_label = f'$\sigma_{{int}})$ {Halpha_label}'
    #
    #         cbar = fig.colorbar(im2)
    #         label_bar = latex_Conver[param]
    #         cbar.set_label(label_bar, rotation=270, labelpad=50, fontsize=20)
    #
    #         ax.update({'title': r'Galaxy {}: {}'.format(obj, title_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    #         ax.set_xlim(100, 210)
    #         ax.set_ylim(90, 240)
    #
    #         # plt.subplots_adjust(top=0.85)
    #         plt.tight_layout()
    #         # plt.show()
    #         plt.savefig(resultsFolder/obj/f'map_{obj}_{dinLabel}_{param}.png')
