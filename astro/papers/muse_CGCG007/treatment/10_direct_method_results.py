import numpy as np
import lime
from lime.plots import STANDARD_PLOT
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import save_log_maps, compute_parameter_distributions, total_abundances_calculation, plot_parameter_image
from lime.io import save_cfg



# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']

fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

# Extensions for fitting files
regions_list = ['MASK_0', 'MASK_1', 'MASK_2']
ionic_param_list = np.array(['n_e', 'T_low', 'cHbeta', 'Ar4', 'Ar3', 'O2', 'O3', 'N2', 'He1', 'S2', 'S3'])
total_abund_list = np.array(['T_high', 'OH', 'NO', 'NH', 'SH', 'ArH', 'S4', 'ICF_S4', 'SO', 'Y_O', 'Y_S', 'S2_S3', 'O2_O3', 'S3_S2',  'eta'])
total_param_list = np.concatenate((ionic_param_list, total_abund_list), axis=0)

# Measurement files and reference files
ref_simulations = ['direct_method']
tech_label = 'neural_fitting'

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = objFolder/'chemistry'
    db_address = objFolder/f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    outputDb = objFolder/f'{obj}_chemical.fits'

    # Header with astronomical data
    hdr = fits.getheader(db_address, f'H1_6563A_flux', ver=1)
    for key, value in hdr.items():
        print(key, value)

    # Loop through the measurement fits
    for j, ref_fit in enumerate(ref_simulations):

        # Generate the map files
        save_log_maps(outputDb, ionic_param_list, chemFolder, maskFits_address, regions_list, ext_log='_CHEMISTRY_OUTPUTS',
                      output_files_prefix=f'{ref_fit}_', page_hdr=hdr)

        # Compute the total abundances and save as fits files
        total_abundances_calculation(ionic_param_list, chemFolder, maskFits_address, regions_list,
                                     ref_fits=f'{ref_fit}_', header=hdr)

        # Plot the parameter value against the image background
        plot_parameter_image(db_address, total_param_list, chemFolder, ref_fit, tech_label=tech_label)

        # Save the parameters
        compute_parameter_distributions(total_param_list, chemFolder, ref_fit, maskFits_address, regions_list, tech_label)


    # # Get mask indeces:
    # spatial_mask_dict = {}
    # with fits.open(maskFits_address) as hdu_masks:
    #     for mask_name in mask_list:
    #         mask_data = hdu_masks[mask_name].data.astype(bool)
    #         spatial_mask_dict[mask_name] = mask_data
    # total_mask = np.array(list(spatial_mask_dict.values()))
    # total_mask = total_mask.sum(axis=0).astype(bool)
    #
    # # # Data for the astronomical coordinates
    # # hdr = fits.getheader(maskFits_address, extname='MASK_0')
    # #
    # # # Generate the map files
    # # save_log_maps(outputDb, param_list, chemFolder, maskFits_address, mask_list, ext_log='_CHEMISTRY_OUTPUTS',
    # #               page_hdr=hdr)
    #
    # # ----------------------------------------- Generate the parameter maps ----------------------------------------
    # flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    # halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    # halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')
    #
    # for param in param_list:
    #
    #     with fits.open(f'{chemFolder}/{param}.fits') as hdu_list:
    #
    #         image_data, image_header = hdu_list[param].data, hdu_list[param].header
    #
    #         defaultConf = STANDARD_PLOT.copy()
    #         rcParams.update(defaultConf)
    #
    #         halpha_cmap = cm.gray.copy()
    #         halpha_cmap.set_under('black')
    #
    #         fig = plt.figure(figsize=(10, 10))
    #         ax = fig.add_subplot(projection=WCS(image_header), slices=('x', 'y'))
    #
    #         bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)
    #
    #         im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
    #         im2 = ax.imshow(image_data)
    #         cbar = fig.colorbar(im2, ax=ax)
    #         param_label = latex_labels[param]
    #         ax.update({'title': r'CGCG007−025, {}'.format(param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    #         ax.set_xlim(120, 210)
    #         ax.set_ylim(110, 220)
    #         plt.savefig(objFolder/f'{obj}_{param_label}_map_directMethod')
    #         # plt.show()
    #
    # # ----------------------------------------- Generate the parameter histograms ----------------------------------------
    # store_dict, err_dict = {}, {}
    # for param in param_list:
    #
    #     with fits.open(f'{chemFolder}/{param}.fits') as hdu_list:
    #
    #         image_data, image_header = hdu_list[param].data, hdu_list[param].header
    #         err_data = hdu_list[f'{param}_err'].data
    #
    #         print(param)
    #
    #         array_data = image_data[total_mask]
    #         array_err = err_data[total_mask]
    #
    #         defaultConf = STANDARD_PLOT.copy()
    #         defaultConf['legend.fontsize'] = 16
    #         rcParams.update(defaultConf)
    #
    #         fig = plt.figure(figsize=(10, 10))
    #         ax = fig.add_subplot()
    #
    #         param_label = latex_labels[param]
    #
    #         # label = r'{} = ${}\pm{}$ {} ({} voxels)'.format(param_label,
    #         #                                              np.round(np.nanmean(array_data), signif_figures[param]),
    #         #                                              np.round(np.nanstd(array_data), signif_figures[param]),
    #         #                                              param_units[param],
    #         #                                              np.sum(total_mask))
    #
    #         median = np.round(np.nanmedian(array_data), signif_figures[param])
    #         upper_limit = np.round(np.nanpercentile(array_data, 84) - np.nanmedian(array_data), signif_figures[param])
    #         lower_limit = np.round(np.nanmedian(array_data) - np.nanpercentile(array_data, 16), signif_figures[param])
    #         n_voxels = np.sum(total_mask)
    #         label = r'{} = ${}^{{{}}}_{{{}}}$ ({} voxels)'.format(param_label, median, upper_limit, lower_limit, n_voxels)
    #         store_dict[f'{param}_array'] = np.array([median, upper_limit, lower_limit, n_voxels])
    #
    #         # Saving pixel error
    #         median_err = np.round(np.nanmedian(array_err), signif_figures[param])
    #         upper_limit_err = np.round(np.nanpercentile(array_err, 84) - np.nanmedian(array_err), signif_figures[param])
    #         lower_limit_err = np.round(np.nanmedian(array_err) - np.nanpercentile(array_err, 16), signif_figures[param])
    #         err_dict[f'{param}_array'] = np.array([median_err, upper_limit_err, lower_limit_err, n_voxels])
    #
    #         ax.hist(array_data, bins=15, label=label)
    #
    #         ax.legend()
    #         ax.update({'title': r'CGCG007−025, {} histogram, direct method'.format(param_label),
    #                    'xlabel': param_label})
    #         plt.savefig(chemFolder/f'{obj}_{param}_histogram_directMethod')
    #         # plt.show()
    #
    #
    # # ----------------------------------------- Generate the O/H, log(N/O) histograms ----------------------------------
    # meta_params = dict.fromkeys(['O2', 'O3', 'N2', 'n_e', 'cHbeta', 'S2', 'S3'], None)
    #
    # for param in meta_params:
    #     with fits.open(f'{chemFolder}/{param}.fits') as hdu_list:
    #         image_data, image_header = hdu_list[param].data, hdu_list[param].header
    #         meta_params[param] = image_data
    #
    # OH_image = 12 + np.log10(np.power(10, meta_params['O2'] - 12) + np.power(10, meta_params['O3'] - 12))
    # OH_array = OH_image[total_mask]
    #
    # NO_image = np.log10(np.power(10, meta_params['N2'] - 12) / np.power(10, meta_params['O2'] - 12))
    # NO_array = NO_image[total_mask]
    #
    # den_ext_image = meta_params['n_e']/meta_params['cHbeta']
    # den_ext_array = den_ext_image[total_mask]
    #
    # O2_O3_image = np.power(10, meta_params['O2'] - 12) / np.power(10, meta_params['O3'] - 12)
    # O2_O3_array = O2_O3_image[total_mask]
    #
    # S2_S3_image = np.power(10, meta_params['S2'] - 12) / np.power(10, meta_params['S3'] - 12)
    # S2_S3_array = S2_S3_image[total_mask]
    #
    # arrays_params = [OH_array, NO_array, den_ext_array, O2_O3_array, S2_S3_array]
    # map_params = [OH_image, NO_image, den_ext_image, O2_O3_image, S2_S3_image]
    #
    # flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    # halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    # halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')
    # for i, param in enumerate(['OH', 'NO', 'nSII_cHbeta', 'O2_O3', 'S2_S3']):
    #
    #     defaultConf = STANDARD_PLOT.copy()
    #     rcParams.update(defaultConf)
    #
    #     halpha_cmap = cm.gray.copy()
    #     halpha_cmap.set_under('black')
    #
    #     fig = plt.figure(figsize=(10, 10))
    #     ax = fig.add_subplot(projection=WCS(image_header), slices=('x', 'y'))
    #
    #     bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)
    #
    #     im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
    #     im2 = ax.imshow(map_params[i])
    #     cbar = fig.colorbar(im2, ax=ax)
    #     param_label = latex_labels[param]
    #     ax.update({'title': r'CGCG007−025, {}, direct method'.format(param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    #     ax.set_xlim(120, 210)
    #     ax.set_ylim(110, 220)
    #     plt.savefig(chemFolder/f'{obj}_{param}_map_directMethod')
    #     # plt.show()
    #
    # for i, param in enumerate(['OH', 'NO', 'nSII_cHbeta', 'O2_O3', 'S2_S3']):
    #
    #     defaultConf = STANDARD_PLOT.copy()
    #     defaultConf['legend.fontsize'] = 16
    #     rcParams.update(defaultConf)
    #
    #     fig = plt.figure(figsize=(10, 10))
    #     ax = fig.add_subplot()
    #     param_label = latex_labels[param]
    #     # label = r'{} = ${}\pm{}$ ({} voxels)'.format(param_label,
    #     #                                              np.round(np.nanmean(arrays_params[i]), signif_figures[param]),
    #     #                                              np.round(np.nanstd(arrays_params[i]), signif_figures[param]),
    #     #                                              np.sum(total_mask))
    #
    #     median = np.round(np.nanmedian(arrays_params[i]), signif_figures[param])
    #     upper_limit = np.round(np.nanpercentile(arrays_params[i], 84) - np.nanmedian(arrays_params[i]), signif_figures[param])
    #     lower_limit = np.round(np.nanmedian(arrays_params[i]) - np.nanpercentile(arrays_params[i], 16), signif_figures[param])
    #     n_voxels = np.sum(total_mask)
    #     label = r'{} = ${}^{{{}}}_{{{}}}$ ({} voxels)'.format(param_label, median, upper_limit, lower_limit, n_voxels)
    #
    #     store_dict[f'{param}_array'] = np.array([median, upper_limit, lower_limit, n_voxels])
    #
    #     ax.hist(arrays_params[i], bins=15, label=label)
    #     ax.legend()
    #     ax.update({'title': r'CGCG007−025, {} histogram, direct method'.format(param_label), 'xlabel': param_label})
    #     plt.savefig(chemFolder/f'{obj}_{param}_histogram_directMethod')
    #     # plt.show()
    #
    #
    # # Save mean values to log
    # save_cfg('../muse_CGCG007.ini', store_dict, section_name='Global_direct_method', clear_section=True)
    # save_cfg('../muse_CGCG007.ini', err_dict, section_name='Global_err_direct_method', clear_section=True)
