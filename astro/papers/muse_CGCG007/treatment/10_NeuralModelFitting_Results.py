import numpy as np
import lime
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, cm, colors, rc_context
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import latex_labels, signif_figures, save_log_maps, convert_dict

from lime.plots import STANDARD_PLOT
from lime.io import save_cfg


def distribution_parametrisation(param_name, param_array):

    median = np.round(np.nanmedian(param_array), signif_figures[param_name])
    up_lim = np.round(np.nanpercentile(param_array, 84) - np.nanmedian(param_array), signif_figures[param_name])
    low_lim = np.round(np.nanmedian(param_array) - np.nanpercentile(param_array, 16), signif_figures[param_name])
    n_voxels = np.sum(~np.isnan(param_array))

    plot_label = r'{} = ${}^{{{}}}_{{{}}}$ ({} voxels)'.format(latex_labels[param_name], median, up_lim, low_lim, n_voxels)
    log_array = np.array([median, up_lim, low_lim, n_voxels])

    return plot_label, log_array


def plot_parameter_image(plot_db_fits, parameter_list, output_folder, conf_label):

    # Image background
    flux6563_image = fits.getdata(plot_db_fits, f'H1_6563A_flux', ver=1)
    halpha_min_level = fits.getval(plot_db_fits, keyword=f'P9050', extname=f'H1_6563A_flux')
    halpha_thresd_level = fits.getval(plot_db_fits, keyword=f'P9250', extname=f'H1_6563A_flux')
    bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)

    # Plot configuration
    defaultConf = STANDARD_PLOT.copy()

    for parameter in parameter_list:

        with fits.open(f'{output_folder}/{parameter}.fits') as hdul:

            param_image, param_hdr = hdul[parameter].data, hdul[parameter].header

            with rc_context(defaultConf):

                halpha_cmap = cm.gray.copy()
                halpha_cmap.set_under('black')

                fig = plt.figure(figsize=(10, 10))
                ax = fig.add_subplot(projection=WCS(param_hdr), slices=('x', 'y'))

                im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
                im2 = ax.imshow(param_image)
                cbar = fig.colorbar(im2, ax=ax)
                plot_wording = {'title': f'CGCG007−025, {latex_labels[parameter]} \n Neural grid sampling, {conf_label}',
                                'xlabel': r'RA', 'ylabel': r'DEC'}
                ax.update(plot_wording)
                ax.set_xlim(120, 210)
                ax.set_ylim(110, 220)
                plt.savefig(chemFolder/f'{parameter}_map_gridSampling_{conf_label}')
                # plt.show()
                plt.close(fig)

    return


def compute_parameter_distributions(parameter_list, output_folder, conf_label, mask_file, mask_list, combined_image=True):

    # Distributions containers
    store_dict, err_dict = {}, {}

    # Plot configuration
    defaultConf = STANDARD_PLOT.copy()
    defaultConf['legend.fontsize'] = 16
    defaultConf['figure.figsize'] = (10, 10)

    # Get regions mask:
    spatial_mask_dict = {}
    with fits.open(mask_file) as hdu_masks:
        for mask_name in mask_list:
            mask_data = hdu_masks[mask_name].data.astype(bool)
            spatial_mask_dict[mask_name] = mask_data

    # region_labels = mask_list[:]
    region_idcs_list = list(spatial_mask_dict.values())

    # Loop throught the parameter file images
    for parameter in parameter_list:

        with fits.open(f'{output_folder}/{conf_label}_{parameter}.fits') as hdu_list:

            # Parameter value and error distribution loop
            for map_type in ['', '_err']:

                param_image = hdu_list[f'{parameter}{map_type}'].data
                dist_container = []
                label_container = []

                # Plot regions distributions
                with rc_context(defaultConf):

                    for i_mask, mask in enumerate(mask_list):

                        data_dist = param_image[region_idcs_list[i_mask]]
                        label_dist, array_dist = distribution_parametrisation(parameter, data_dist)

                        label_container.append(label_dist)
                        dist_container.append(data_dist)

                        # Store the distribution
                        ref_dist = f'{convert_dict[parameter]}{map_type}_{mask}'
                        store_dict[ref_dist] = array_dist

                    # Make the plot
                    fig, ax = plt.subplots()
                    ax.hist(dist_container, bins=15, label=label_container, stacked=True)
                    ax.set_title(f'CGCG007−025, {latex_labels[parameter]} \n Neural model fitting {conf_label} \n regions histogram')
                    ax.set_xlabel(latex_labels[parameter])
                    ax.legend()
                    plt.savefig(chemFolder/f'neuralFitting_{conf_label}_{parameter}_regions_histogram{map_type}')
                    plt.close(fig)
                    # plt.show()

                # Plot global distribution
                total_idcs_mask = np.array(region_idcs_list).sum(axis=0).astype(bool)
                ref_global = 'global'

                # Plot total distribution
                with rc_context(defaultConf):

                    data_dist = param_image[total_idcs_mask]
                    label_dist, array_dist = distribution_parametrisation(parameter, param_image)

                    # Store the distribution
                    ref_dist = f'{convert_dict[parameter]}{map_type}_{ref_global}'
                    store_dict[ref_dist] = array_dist

                    # Make the plot
                    fig, ax = plt.subplots()
                    ax.hist(data_dist, bins=15, label=label_dist)
                    ax.set_title(f'CGCG007−025, {latex_labels[parameter]} \n Neural model fitting {conf_label} \n all voxels histogram')
                    ax.set_xlabel(latex_labels[parameter])
                    ax.legend()
                    plt.savefig(chemFolder / f'neuralFitting_{conf_label}_{parameter}_global_histogram{map_type}')
                    plt.close(fig)

    # Save to a file
    save_cfg('../muse_CGCG007.ini', store_dict, section_name=f'GridSampling_{conf_label}', clear_section=True)

    return


# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']

fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

# Extensions for fitting files
regions_list = ['MASK_0', 'MASK_1', 'MASK_2']
param_list = np.array(['logOH', 'logNO', 'logU'])

# Measurement files and reference files
ref_simulations = ['localErr', 'HIICHImistry', 'noOII', 'maxErr']

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = objFolder/'chemistry'
    db_address = objFolder/f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Header with astronomical data
    hdr = fits.getheader(maskFits_address, extname='MASK_0')

    # Loop through the measurement fits
    for j, ref_fit in enumerate(ref_simulations):

        # # Generate the parameter maps
        # grid_fits_file = objFolder/f'{obj}_grid_sampling_{ref_fit}.fits'
        # save_log_maps(grid_fits_file, param_list, chemFolder, maskFits_address, regions_list,
        #               output_files_prefix=f'{ref_fit}_', ext_log='_GRIDSAMPLER_OUTPUTS', page_hdr=hdr)

        # Plot the parameter value against the image background
        # plot_parameter_image(db_address, param_list, chemFolder, ref_fit)

        # Generate the parameter maps and store the distributions
        compute_parameter_distributions(param_list, chemFolder, ref_fit, maskFits_address, regions_list)
