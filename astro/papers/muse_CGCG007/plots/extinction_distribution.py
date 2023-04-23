import numpy as np
import lime

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import chemical_lines_indexing, import_muse_fits
from astropy.io import fits
from astropy.table import Table
from matplotlib import pyplot as plt, rcParams, cm, colors, rc_context
from astropy.wcs import WCS
from src.specsiser.plots import latex_labels


# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
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

merge_dict = {'O2_7319A_b': 'O2_7319A-O2_7330A'}

R_v = obsData['Extinction']['R_v']
red_law = obsData['Extinction']['red_law']
image_size = obsData['sample_data']['grid_shape_array'].astype(int)

lines_mask_dict = {'MASK_0': ['H1_4861A', 'H1_9229A', 'H1_9015A', 'H1_8863A', 'H1_8750A'],
                   'MASK_1': ['H1_4861A', 'H1_9229A', 'H1_9015A', 'H1_8863A'],
                   'MASK_2': ['H1_4861A', 'H1_6563A', 'H1_9229A', 'H1_9015A'],
                   'MASK_3': ['H1_4861A', 'H1_6563A'],
                   'MASK_4': ['H1_4861A', 'H1_6563A'],
                   'MASK_5': ['H1_4861A', 'H1_6563A']}



STANDARD_PLOT = {'figure.figsize': (12, 12),
                 'axes.titlesize': 18,
                 'axes.labelsize': 30,
                 'legend.fontsize': 18,
                 'xtick.labelsize': 25,
                 'ytick.labelsize': 25}

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = f'{objFolder}/chemistry'
    obsLog_addresss = objFolder / f'{obj}_linesLog.fits'
    absLog_addresss = objFolder / f'{obj}_linesLog_abs.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    db_address = objFolder / f'{obj}_database.fits'

    # Output files
    output_map = objFolder/f'{obj}_HI_extinction.fits'
    images_dict = {'cHbeta': np.full(image_size, np.nan),
                   'cHbeta_err': np.full(image_size, np.nan)}

    # # ----------------------------------------- Generate the image plots ----------------------------------------
    flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')

    param = 'cHbeta'
    with fits.open(output_map) as hdu_list:
        image_data, image_header = hdu_list[param].data, hdu_list[param].header

    mask_sum = np.full(image_size, False)
    array_container = []
    region_list = [0, 1, 2, 3, 4, 5]
    voxel_count = np.zeros(len(region_list))
    cHbeta_array = [None] * len(voxel_count)
    for idx_region in region_list:

        # Voxel mask
        region_label = f'mask_{idx_region}'
        region_mask = fits.getdata(maskFits_address, region_label, ver=1)
        region_mask = region_mask.astype(bool)
        array_container.append(image_data[region_mask])
        mask_sum += region_mask
        voxel_count[idx_region] = np.sum(region_mask)
        cHbeta_array[idx_region] = (np.nanmean(image_data[region_mask]), np.nanstd(image_data[region_mask]))


    with rc_context(STANDARD_PLOT):

        # halpha_cmap = cm.gray
        # halpha_cmap.set_under('black')

        # fig = plt.figure(figsize=(10, 10))
        # ax = fig.add_subplot(projection=WCS(image_header), slices=('x', 'y'))
        #
        # bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)
        # im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
        #
        # param_min, param_max = np.nanmin(image_data), np.nanmax(image_data)
        # divnorm = colors.TwoSlopeNorm(vcenter=0.0, vmin=-0.05, vmax=0.5, )
        #
        # # im2 = ax.imshow(image_data, norm=colors.LogNorm(vmin=0.01, vmax=0.8))
        # image_data[~mask_sum] = np.nan
        # im2 = ax.imshow(image_data, norm=divnorm)
        #
        # cbar = fig.colorbar(im2, ax=ax)
        # param_label = latex_labels[param]
        # ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
        # ax.set_xlim(95, 210)
        # ax.set_ylim(80, 222)
        # plt.show()
        # plt.savefig(objFolder/f'{obj}_extinction_map.png')

        fig = plt.figure(dpi=600)
        ax = fig.add_subplot()

        # ax.hist(image_data[mask_sum], bins=50, log=True)
        cmap = cm.get_cmap(name='viridis_r')
        colorNorm = colors.Normalize(0, len(region_list))
        colors_hist = [cmap(colorNorm(i)) for i in region_list]

        cHbeta_label = [r'$cH(\beta) = {:0.2f}\pm{:0.2f}$'.format(cHbeta_array[i][0], cHbeta_array[i][1]) for i in region_list]
        labels_hist = [f'{cHbeta_label[i]} (region {i}, {voxel_count[i]:.0f} voxels)' for i in region_list]
        ax.hist(array_container, label=labels_hist, bins=50, log=True, stacked=True, color=colors_hist)

        ax.update({'xlabel': r'$cH(\beta)$', 'ylabel': r'Voxel count'})
        # ax.legend()
        # plt.show()
        plt.savefig(plotsFolder/'CGCG007_extinction_distribution.png')