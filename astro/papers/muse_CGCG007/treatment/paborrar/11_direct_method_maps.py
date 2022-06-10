import numpy as np
import lime
from lime.plots import STANDARD_PLOT
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import save_log_maps, latex_labels, signif_figures, param_units
from lime.io import save_cfg

# Declare data and files location
obsData = lime.load_cfg('../../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']

fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

voxel_grid_size = obsData['sample_data']['grid_shape_array']

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    maskFits_address, mask_list = objFolder/f'{obj}_masks.fits', ['MASK_0', 'MASK_1', 'MASK_2']
    db_address = objFolder/f'{obj}_database.fits'
    outputDb = objFolder/f'{obj}_chemical.fits'
    chemFolder = objFolder/'chemistry'

    # Parameters to plot
    param_list = np.array(['n_e', 'T_low', 'cHbeta', 'Ar4', 'Ar3', 'O2', 'O3', 'N2', 'He1', 'S2', 'S3'])

    # Get mask indeces:
    spatial_mask_dict = {}
    with fits.open(maskFits_address) as hdu_masks:
        for mask_name in mask_list:
            spatial_mask_dict[mask_name] = hdu_masks[mask_name].data.astype(bool)
    total_mask = np.array(list(spatial_mask_dict.values()))
    total_mask = total_mask.sum(axis=0).astype(bool)

    # # Data for the astronomical coordinates
    # hdr = fits.getheader(maskFits_address, extname='MASK_0')
    #
    # # Generate the map files
    # save_log_maps(outputDb, param_list, chemFolder, maskFits_address, mask_list, ext_log='_CHEMISTRY_OUTPUTS',
    #               page_hdr=hdr)

    # ----------------------------------------- Generate the parameter histograms ----------------------------------------
    for param in param_list:

        with fits.open(f'{chemFolder}/{param}.fits') as hdu_list:

            image_data, image_header = hdu_list[param].data, hdu_list[param].header

            defaultConf = STANDARD_PLOT.copy()
            defaultConf['legend.fontsize'] = 16
            rcParams.update(defaultConf)

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot()

            param_label = latex_labels[param]

            store_dict = {}
            array_container, data_labels = [], []
            for i_mask, mask_items in enumerate(spatial_mask_dict.items()):
                mask_name, mask_image = mask_items
                array_data = image_data[mask_image]

                median = np.round(np.nanmedian(array_data), signif_figures[param])
                upper_limit = np.round(np.nanpercentile(array_data, 84) - np.nanmedian(array_data), signif_figures[param])
                lower_limit = np.round(np.nanmedian(array_data) - np.nanpercentile(array_data, 16), signif_figures[param])
                n_voxels = np.sum(mask_image)
                label = r'{} = ${}^{{{}}}_{{{}}}$ ({} voxels)'.format(param_label, median, upper_limit, lower_limit, n_voxels)

                array_container.append(array_data)
                data_labels.append(label)

                store_dict = {f'{param}_array': np.array([median, upper_limit, lower_limit, n_voxels])}
                save_cfg('../../muse_CGCG007.ini', store_dict, section_name=f'{mask_name}_direct_method')

            ax.hist(array_container, bins=15, label=data_labels, stacked=True)

            ax.legend()
            ax.update({'title': r'CGCG007−025, {} histogram, direct method'.format(param_label), 'xlabel': param_label})
            plt.show()
            # plt.savefig(chemFolder/f'{obj}_{param}_histogram_directMethod')

    # ----------------------------------------- Generate the O/H, log(N/O) histograms ----------------------------------
    meta_params = dict.fromkeys(['O2', 'O3', 'N2', 'n_e', 'cHbeta', 'S2', 'S3'], None)

    for param in meta_params:
        with fits.open(f'{chemFolder}/{param}.fits') as hdu_list:
            image_data, image_header = hdu_list[param].data, hdu_list[param].header
            meta_params[param] = image_data

    OH_image = 12 + np.log10(np.power(10, meta_params['O2'] - 12) + np.power(10, meta_params['O3'] - 12))
    OH_array = OH_image[total_mask]

    NO_image = np.log10(np.power(10, meta_params['N2'] - 12) / np.power(10, meta_params['O2'] - 12))
    NO_array = NO_image[total_mask]

    den_ext_image = meta_params['n_e']/meta_params['cHbeta']
    den_ext_array = den_ext_image[total_mask]

    O2_O3_image = np.power(10, meta_params['O2'] - 12) / np.power(10, meta_params['O3'] - 12)
    O2_O3_array = O2_O3_image[total_mask]

    S2_S3_image = np.power(10, meta_params['S2'] - 12) / np.power(10, meta_params['S3'] - 12)
    S2_S3_array = S2_S3_image[total_mask]

    map_params = [OH_image, NO_image, den_ext_image, O2_O3_image, S2_S3_image]

    for i, param in enumerate(['OH', 'NO', 'nSII_cHbeta', 'O2_O3', 'S2_S3']):

        defaultConf = STANDARD_PLOT.copy()
        defaultConf['legend.fontsize'] = 16
        rcParams.update(defaultConf)

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot()
        param_label = latex_labels[param]

        array_container, data_labels = [], []
        for i_mask, mask_items in enumerate(spatial_mask_dict.items()):
            mask_name, mask_image = mask_items
            array_data = map_params[i][mask_image]

            median = np.round(np.nanmedian(array_data), signif_figures[param])
            upper_limit = np.round(np.nanpercentile(array_data, 84) - np.nanmedian(array_data), signif_figures[param])
            lower_limit = np.round(np.nanmedian(array_data) - np.nanpercentile(array_data, 16), signif_figures[param])
            n_voxels = np.sum(mask_image)
            label = r'{} = ${}^{{{}}}_{{{}}}$ ({} voxels)'.format(param_label, median, upper_limit, lower_limit, n_voxels)

            array_container.append(array_data)
            data_labels.append(label)

            store_dict = {f'{param}_array': np.array([median, upper_limit, lower_limit, n_voxels])}
            save_cfg('../../muse_CGCG007.ini', store_dict, section_name=f'{mask_name}_direct_method')

        ax.hist(array_container, bins=15, label=data_labels, stacked=True)
        ax.legend()
        ax.update({'title': r'CGCG007−025, {} histogram, direct method'.format(param_label), 'xlabel': param_label})
        # plt.savefig(chemFolder/f'{obj}_{param}_histogram_directMethod')
        plt.show()


