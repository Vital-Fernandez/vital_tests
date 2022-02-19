import src.specsiser as sr
from delete.data_printing import latex_labels
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from src.specsiser.tools.line_measure import STANDARD_PLOT


# Declare data and files location
obsData = sr.loadConfData('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

voxel_grid_size = obsData['sample_data']['grid_shape_array']

verbose = True

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_address = objFolder / f'{obj}_database.fits'
    fits_results_file = resultsFolder/obj/f'{obj}_chemical.fits'

    # Output files
    chemMaps_fits_address = objFolder/f'{obj}_ChemParamMaps.fits'
    theoLineMaps_fits_address = objFolder/f'{obj}_fitTheoFluxesMaps.fits'

    # # ----------------------------------------- Generate the image data -----------------------------------------
    # 
    # # Empty containers for the images
    # image_dict = {}
    # theoLines_dict = {}
    # 
    # # Open the chemical results file
    # with fits.open(fits_results_file) as hdul:
    # 
    #     # Loop throught the line regions
    #     for idx_region in [0, 1, 2]:
    # 
    #         region_label = f'region_{idx_region}'
    #         region_mask = fits.getdata(db_address, region_label, ver=1)
    #         region_mask = region_mask.astype(bool)
    #         idcs_voxels = np.argwhere(region_mask)
    # 
    #         print(f'- Importing data {region_label}')
    # 
    #         # Loop through the region voxels
    #         for idx_voxel, idx_pair in enumerate(idcs_voxels):
    #             idx_j, idx_i = idx_pair
    # 
    #             outputsLabel = f'{idx_j}-{idx_i}_CHEMISTRY_OUTPUTS'
    #             tracesLabel = f'{idx_j}-{idx_i}_CHEMISTRY_TRACES'
    # 
    #             # Load lines log data and store it as an image
    #             if outputsLabel in hdul:
    #                 chem_traces = hdul[tracesLabel].data
    #                 chem_labels = hdul[outputsLabel].data['parameters_list']
    # 
    #                 for param in chem_labels:
    #                     trace_i = chem_traces[param]
    #                     trace_mean, trace_std = np.mean(trace_i), np.std(trace_i)
    # 
    #                     # Create image map if not available
    #                     if param not in image_dict:
    #                         image_dict[param] = np.full(voxel_grid_size.astype(int), np.nan)
    #                         image_dict[f'{param}_err'] = np.full(voxel_grid_size.astype(int), np.nan)
    # 
    #                     # Chemical parameter case
    #                     image_dict[param][idx_j, idx_i] = trace_mean
    #                     image_dict[f'{param}_err'][idx_j, idx_i] = trace_std
    # 
    # 
    # # Storing the dictionaries as a fits image files
    # dict_list = [image_dict]
    # file_list = [chemMaps_fits_address]
    # 
    # for j, data_dict in enumerate(dict_list):
    #     new_hdul = fits.HDUList()
    #     new_hdul.append(fits.PrimaryHDU())
    # 
    #     # Second page for the fits file plot configuration
    #     hdr_plot = fits.getheader(db_address, extname='PlotConf')
    #     hdu_table = fits.BinTableHDU.from_columns(columns=[], header=hdr_plot, name='PlotConf')
    #     new_hdul.append(hdu_table)
    # 
    #     for param_line, param_map in data_dict.items():
    #         print(param_line)
    #         new_hdul.append(fits.ImageHDU(name=param_line, data=param_map, ver=1, header=fits.Header({'param': param_line})))
    #     new_hdul.writeto(file_list[j], overwrite=True, output_verify='fix')

    # ----------------------------------------- Generate the image plots ----------------------------------------
    hdr_plot = fits.getheader(chemMaps_fits_address, extname='PlotConf')
    flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')


    with fits.open(chemMaps_fits_address) as hdu_list:

        for ext in hdu_list:

            ext_label = ext.name

            if (ext_label not in ['PRIMARY', 'PLOTCONF']) and ('_ERR' not in ext_label):

                image_data = hdu_list[ext_label].data
                image_header = hdu_list[ext_label].header

                defaultConf = STANDARD_PLOT.copy()
                rcParams.update(defaultConf)

                halpha_cmap = cm.gray
                halpha_cmap.set_under('black')

                fig = plt.figure(figsize=(10, 10))
                ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))

                param_label = image_header['param']
                # param_image = fits.getdata(chemMaps_fits_address, param_label, ver=1)

                im = ax.imshow(flux6563_image, cmap=halpha_cmap,
                               norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))
                im2 = ax.imshow(image_data)

                cbar = fig.colorbar(im2, ax=ax)

                param_label = latex_labels[param_label]
                ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
                ax.set_xlim(120, 210)
                ax.set_ylim(110, 220)
                plt.tight_layout()
                # plt.savefig(objFolder/f'{obj}_parameter_map_{ext_label}')
                plt.show()
