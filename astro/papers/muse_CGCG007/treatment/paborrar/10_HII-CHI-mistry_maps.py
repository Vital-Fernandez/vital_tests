import numpy as np
import lime
from lime.plots import STANDARD_PLOT
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import latex_labels, signif_figures
from lime.io import save_cfg

# Declare data and files location
obsData = lime.load_cfg('../../muse_CGCG007.ini')
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

    # Parameters to plot
    param_list = np.array(['OH', 'NO', 'logU'])

    # Get mask indeces:
    spatial_mask_dict = {}
    with fits.open(maskFits_address) as hdu_masks:
        for mask_name in mask_list:
            mask_data = hdu_masks[mask_name].data.astype(bool)
            spatial_mask_dict[mask_name] = mask_data
    total_mask = np.array(list(spatial_mask_dict.values()))
    total_mask = total_mask.sum(axis=0).astype(bool)

    # ----------------------------------------- Generate the parameter maps ----------------------------------------
    flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')

    for param in param_list:

        with fits.open(f'{chemFolder}/{param}_HII-CHI-mistry_output.fits') as hdu_list:

            image_data, image_header = hdu_list[param].data, hdu_list[param].header

            idcs_data = ~np.isnan(image_data)

            defaultConf = STANDARD_PLOT.copy()
            rcParams.update(defaultConf)

            halpha_cmap = cm.gray.copy()
            halpha_cmap.set_under('black')

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(projection=WCS(image_header), slices=('x', 'y'))

            bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)

            im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
            im2 = ax.imshow(image_data)
            cbar = fig.colorbar(im2, ax=ax)
            param_label = latex_labels[param]
            ax.update({'title': r'CGCG007−025, {}, HII-CHI-mistry'.format(param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
            ax.set_xlim(120, 210)
            ax.set_ylim(110, 220)
            plt.savefig(chemFolder/f'{obj}_{param}_map_HII-CHI-mistry')
            # plt.show()

    # ----------------------------------------- Generate the parameter histograms ----------------------------------------
    store_dict, err_dict = {}, {}
    for param in param_list:

        with fits.open(f'{chemFolder}/{param}_HII-CHI-mistry_output.fits') as hdu_list:

            image_data, image_header, image_error = hdu_list[param].data, hdu_list[param].header, hdu_list[f'{param}_err'].data
            array_data = image_data[total_mask]
            array_err = image_error[total_mask]

            defaultConf = STANDARD_PLOT.copy()
            defaultConf['legend.fontsize'] = 16
            rcParams.update(defaultConf)

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot()

            param_label = latex_labels[param]

            # label = r'{} = ${}\pm{}$ {} ({} voxels)'.format(param_label,
            #                                              np.round(np.nanmean(array_data), signif_figures[param]),
            #                                              np.round(np.nanstd(array_data), signif_figures[param]),
            #                                              param_units[param],
            #                                              np.sum(total_mask))

            # Measurement distribution
            median = np.round(np.nanmedian(array_data), signif_figures[param])
            upper_limit = np.round(np.nanpercentile(array_data, 84) - np.nanmedian(array_data), signif_figures[param])
            lower_limit = np.round(np.nanmedian(array_data) - np.nanpercentile(array_data, 16), signif_figures[param])
            n_voxels = np.sum(total_mask)
            label = r'{} = ${}^{{{}}}_{{{}}}$ ({} voxels)'.format(param_label, median, upper_limit, lower_limit, n_voxels)

            # Error distribution
            median_err = np.round(np.nanmedian(array_err), signif_figures[param])
            upper_limit_err = np.round(np.nanpercentile(array_err, 84) - np.nanmedian(array_err), signif_figures[param])
            lower_limit_err = np.round(np.nanmedian(array_err) - np.nanpercentile(array_err, 16), signif_figures[param])
            n_voxels_err = np.sum(total_mask)

            store_dict[f'{param}_array'] = np.array([median, upper_limit, lower_limit, n_voxels])
            err_dict[f'{param}_array'] = np.array([median_err, upper_limit_err, lower_limit_err, n_voxels_err])

            ax.hist(array_data, bins=15, label=label)

            ax.legend()
            ax.update({'title': r'CGCG007−025, {} histogram, HII-CHI-mistry'.format(param_label),
                       'xlabel': param_label})
            plt.savefig(chemFolder/f'{obj}_{param}_histogram_HII-CHI-mistry')
            # plt.show()

    # Save mean values to log
    save_cfg('../../muse_CGCG007.ini', store_dict, section_name='Global_HII-CHI-mistry', clear_section=True)
    save_cfg('../../muse_CGCG007.ini', err_dict, section_name='Global_err_HII-CHI-mistry', clear_section=True)
