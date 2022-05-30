import numpy as np
import lime
from lime.plots import STANDARD_PLOT
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import latex_labels, signif_figures, save_log_maps, convert_dict
from lime.io import save_cfg

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']

fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    maskFits_address, mask_list = objFolder/f'{obj}_masks.fits', ['MASK_0', 'MASK_1', 'MASK_2']
    db_address = objFolder/f'{obj}_database.fits'
    chemFolder = objFolder/'chemistry'
    grid_fits_file = objFolder/f'{obj}_grid_sampling_maxErr.fits'

    # Parameters to plot
    param_list = np.array(['logOH', 'logNO', 'logU'])

    # Get mask indeces:
    spatial_mask_dict = {}
    with fits.open(maskFits_address) as hdu_masks:
        for mask_name in mask_list:
            mask_data = hdu_masks[mask_name].data.astype(bool)
            spatial_mask_dict[mask_name] = mask_data
    total_mask = np.array(list(spatial_mask_dict.values()))
    total_mask = total_mask.sum(axis=0).astype(bool)

    # # Data for the astronomical coordinates
    # hdr = fits.getheader(maskFits_address, extname='MASK_0')
    #
    # # Generate the map files
    # save_log_maps(grid_fits_file, param_list, chemFolder, maskFits_address, mask_list, ext_log='_GRIDSAMPLER_OUTPUTS',
    #               page_hdr=hdr)

    # ----------------------------------------- Generate the parameter histograms ----------------------------------------

    for param in param_list:

        with fits.open(f'{chemFolder}/MAX_ERR{param}.fits') as hdu_list:

            image_data, image_header = hdu_list[f'{param}_err'].data, hdu_list[f'{param}_err'].header
            # array_data = image_data[total_mask]

            defaultConf = STANDARD_PLOT.copy()
            defaultConf['legend.fontsize'] = 16
            rcParams.update(defaultConf)

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot()

            param_label = latex_labels[param]

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

                store_dict = {f'{convert_dict[param]}_array': np.array([median, upper_limit, lower_limit, n_voxels])}
                save_cfg('../muse_CGCG007.ini', store_dict, section_name=f'{mask_name}_err_GridSampling_MaxErr')

            ax.hist(array_container, bins=15, label=data_labels, stacked=True)

            ax.legend()
            ax.update({'title': r'CGCG007−025, {} uncertainty histogram, Grid sampling, Max Err'.format(param_label), 'xlabel': param_label})
            plt.savefig(chemFolder/f'{obj}_{param}_err_histogram_GridSamplingMAXERR')
            # plt.show()

