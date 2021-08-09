import numpy as np
import src.specsiser as sr
from src.specsiser.data_printing import latex_labels
from pathlib import Path
from astro.data.muse.common_methods import background_color, DARK_PLOT, label_Conver, latex_Conver, dinamicLines
from timeit import default_timer as timer
from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS


# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

voxel_grid_size = obsData['sample_data']['grid_shape_array']

verbose = True

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        objFolder = resultsFolder/obj
        voxelFolder = resultsFolder/obj/'voxel_data'
        db_address = objFolder / f'{obj}_database.fits'
        chem_conf_file = dataFolder/obj/'chemical_model_config.txt'
        fits_results_file = resultsFolder/obj/f'{obj}_chemical.fits'

        # Load chemical conf
        chem_conf = sr.loadConfData(chem_conf_file, group_variables=False)

        # Output files
        chemMaps_fits_address = objFolder/f'{obj}_ChemParamMaps.fits'
        theoLineMaps_fits_address = objFolder/f'{obj}_fitTheoFluxesMaps.fits'

        # # ----------------------------------------- Generate the image data -----------------------------------------
        #
        # # Empty containers for the images
        # image_dict = {}
        # theoLines_dict = {}
        # chem_params = chem_conf['inference_model_configuration']['parameter_list']
        # for param in chem_params:
        #     image_dict[param] = np.full(voxel_grid_size.astype(int), np.nan)
        #     image_dict[f'{param}_err'] = np.full(voxel_grid_size.astype(int), np.nan)
        #
        # # Open the chemical results file
        # with fits.open(fits_results_file) as hdul:
        #
        #     # Loop throught the line regions
        #     for idx_region in [0, 1, 2]:
        #         region_label = f'region_{idx_region}'
        #         region_mask = fits.getdata(db_address, region_label, ver=1)
        #         region_mask = region_mask.astype(bool)
        #         idcs_voxels = np.argwhere(region_mask)
        #         print(f'- Importing data {region_label}')
        #
        #         # Loop through the region voxels
        #         for idx_voxel, idx_pair in enumerate(idcs_voxels):
        #             idx_j, idx_i = idx_pair
        #             logLabel = f'{idx_j}-{idx_i}_CHEMISTRY'
        #
        #             # Load lines log data and store it as an image
        #             if logLabel in hdul:
        #                 chem_traces = hdul[logLabel].data
        #                 chem_labels = chem_traces.names
        #
        #                 for param in chem_labels:
        #                     trace_i = chem_traces[param]
        #                     trace_mean, trace_std = np.mean(trace_i), np.std(trace_i)
        #
        #                     # Chemical parameter case
        #                     if param in chem_params:
        #                         image_dict[param][idx_j, idx_i] = trace_mean
        #                         image_dict[f'{param}_err'][idx_j, idx_i] = trace_std
        #
        #                     # Emission line case
        #                     else:
        #                         if param not in image_dict:
        #                             theoLines_dict[f'{param}'] = np.full(voxel_grid_size.astype(int), np.nan)
        #                             theoLines_dict[f'{param}_err'] = np.full(voxel_grid_size.astype(int), np.nan)
        #
        #                         theoLines_dict[param][idx_j, idx_i] = trace_mean
        #                         theoLines_dict[f'{param}_err'][idx_j, idx_i] = trace_std
        #
        # # Storing the dictionaries as a fits image files
        # dict_list = [image_dict, theoLines_dict]
        # file_list = [chemMaps_fits_address, theoLineMaps_fits_address]
        #
        # for i, data_dict in enumerate(dict_list):
        #     new_hdul = fits.HDUList()
        #     new_hdul.append(fits.PrimaryHDU())
        #
        #     # Second page for the fits file plot configuration
        #     hdr_plot = fits.getheader(db_address, extname='PlotConf')
        #     hdu_table = fits.BinTableHDU.from_columns(columns=[], header=hdr_plot, name='PlotConf')
        #     new_hdul.append(hdu_table)
        #
        #     for param_line, param_map in data_dict.items():
        #         new_hdul.append(fits.ImageHDU(name=param_line, data=param_map, ver=1))
        #     new_hdul.writeto(file_list[i], overwrite=True, output_verify='fix')


        # ----------------------------------------- Generate the image plots ----------------------------------------
        hdr_plot = fits.getheader(chemMaps_fits_address, extname='PlotConf')
        flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
        halpha_min_level = fits.getval(db_address, keyword=f'P9000', extname=f'H1_6563A_flux')
        halpha_thresd_level = fits.getval(db_address, keyword=f'P8000', extname=f'H1_6563A_flux')

        defaultConf = DARK_PLOT.copy()
        rcParams.update(defaultConf)

        halpha_cmap = cm.gray
        halpha_cmap.set_under(background_color)

        # param lines
        for i, param in enumerate(['T_low', 'n_e', 'cHbeta', 'N/O', 'O']):

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))

            if param != 'N/O' and param  != 'O':
                param_image = fits.getdata(chemMaps_fits_address, param, ver=1)
            if param == 'O':
                O2 = np.power(10, fits.getdata(chemMaps_fits_address, 'O2', ver=1)-12)
                O3 = np.power(10, fits.getdata(chemMaps_fits_address, 'O3', ver=1)-12)
                param_image = 12 + np.log10(O2 + O3)
                latex_labels['O'] = r'12+log(O)'

            else:
                N2 = np.power(10, fits.getdata(chemMaps_fits_address, 'N2', ver=1)-12)
                O2 = np.power(10, fits.getdata(chemMaps_fits_address, 'O2', ver=1)-12)
                param_image = np.log10(N2/O2)
                latex_labels['N/O'] = r'log(N/O)'
            im = ax.imshow(flux6563_image, cmap=halpha_cmap,
                           norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))
            im2 = ax.imshow(param_image)

            cbar = fig.colorbar(im2, ax=ax)

            param_label = latex_labels[param]
            ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
            plt.tight_layout()
            plt.show()



        # # Loop throught the line regions
        # for idx_region in [0, 1, 2]:
        #
        #     # Voxel mask
        #     region_label = f'region_{idx_region}'
        #     region_mask = fits.getdata(db_addresss, region_label, ver=1)
        #     region_mask = region_mask.astype(bool)
        #     idcs_voxels = np.argwhere(region_mask)
        #     n_voxels = idcs_voxels.shape[0]
        #     print(f'\nTreating {region_label} consisting of {n_voxels}. The estimated time is {(n_voxels*1.65)/60:0.1f} hrs')
        #
        #     for idx_voxel, idx_pair in enumerate(idcs_voxels):
        #
        #         idx_j, idx_i = idx_pair
        #
        #
        # # Declare voxels to analyse
        # flux6312_image = fits.getdata(db_addresss, f'{ref_flux_line}_flux', ver=1)
        # flux6312_levels = np.nanpercentile(flux6312_image, pertil_array)
        #
        # flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
        # flux6563_levels = np.nanpercentile(flux6563_image, pertil_array)
        #
        # # Search within that limit
        # maFlux_image = np.ma.masked_where((flux6312_image >= flux6312_levels[sulfur_bdry]) &
        #                                   (flux6563_image > flux6563_levels[hydrogen_bdry]),
        #                                   flux6563_image)
        # idcs_voxels = np.argwhere(maFlux_image.mask)
        # n_voxels = idcs_voxels.shape[0]
        #
        # if verbose:
        #     fig = plt.figure(figsize=(12, 8))
        #     ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
        #     ax.update({'title': r'{} galaxy, $H\alpha$ flux'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
        #     im = ax.imshow(maFlux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=flux6563_levels[1],
        #                    vmin=flux6563_levels[1], base=10))
        #     cntr1 = ax.contour(flux6312_image, levels=[flux6312_levels[sulfur_bdry]], colors='yellow', alpha=0.5)
        #     cntr2 = ax.contour(flux6563_image, levels=[flux6563_levels[hydrogen_bdry]], colors='red', alpha=0.5)
        #     plt.show()
        #
        # print(f'\nUsing line [SIII]6312 at percentile {pertil_array[sulfur_bdry]} = {flux6312_levels[sulfur_bdry]:.2f}'
        #       f' ({idcs_voxels.shape[0]} pixels)')
        #
        # parameter_list = ['n_e', 'T_low', 'S3', 'cHbeta']
        # parameter_image = np.empty(flux6563_image.shape)
        #
        # # Construct the images
        # maps_dict = {}
        # voxel_n = 0
        # with fits.open(fits_model_file) as hdul:
        #
        #     for parameter in parameter_list:
        #
        #         parameter_image = np.empty(flux6563_image.shape)
        #         parameter_image[:] = np.nan
        #
        #         for idx_voxel, idx_pair in enumerate(idcs_voxels):
        #
        #             idx_j, idx_i = idx_pair
        #
        #             chem_ref = f'{idx_j}-{idx_i}_chemistry'
        #
        #             if chem_ref in hdul:
        #                 header = hdul[chem_ref].header
        #                 data = hdul[chem_ref].data
        #                 trace = data[parameter]
        #                 parameter_image[idx_j, idx_i] = trace.mean()
        #                 voxel_n +=1
        #
        #         maps_dict[parameter] = parameter_image
        #
        # # Plot the maps
        # print(f'Number of voxels treated: {voxel_n/4}/{n_voxels}')
        # for parameter, param_map in maps_dict.items():
        #
        #     fig = plt.figure(figsize=(12, 8))
        #     ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
        #     im = ax.imshow(maFlux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=flux6563_levels[1],
        #                                                                       vmin=flux6563_levels[1], base=10))
        #     im2 = ax.imshow(param_map)
        #     cbar = fig.colorbar(im2)
        #     ax.update({'title': r'{} galaxy: {}'.format(obj, latex_labels[parameter]), 'xlabel': r'RA', 'ylabel': r'DEC'})
        #     plt.tight_layout()
        #     # plt.savefig(resultsFolder/obj/f'map_{obj}_{parameter}.png', )
        #     plt.show()
