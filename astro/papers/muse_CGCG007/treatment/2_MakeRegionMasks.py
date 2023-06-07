import numpy as np
import pandas as pd
import lime as lm
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, cm, colors, patches, rcParams
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits
from lime.plots import STANDARD_PLOT

# Declare data and files location
obsData = lm.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

dict_errs = {}
dict_nan_values = {}

ref_flux_line = 'S3_6312A'

for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_addresss = objFolder/f'{obj}_database.fits'
    mask_address = dataFolder/f'{obj}_mask.txt'

    # Output data:
    hdul_lineslog = fits.HDUList()
    masks_fits = objFolder/f'{obj}_masks.fits'
    masks_plot = objFolder/f'{obj}_masks.png'

    # Load data
    wave, cube, header = import_muse_fits(cube_address)
    wave_rest = wave / (1 + z_objs[i])
    mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)

    # Declare voxels to analyse
    invers_pertil_array = pertil_array[::-1]
    flux6312_image = fits.getdata(db_addresss, f'{ref_flux_line}_flux', ver=1)
    flux6312_levels = np.nanpercentile(flux6312_image, invers_pertil_array)

    flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
    flux6563_levels = np.nanpercentile(flux6563_image, invers_pertil_array)

    Halpha_idxMinPercentil = 6
    Halpha_min_level = flux6563_levels[5]
    SIII_idxMinPercentil = 3

    for i in np.arange(len(flux6312_levels)):
        print(f'{i}) S3_6312A band image {invers_pertil_array[i]} percentil : flux {flux6312_levels[i]:.1f} for {np.sum(flux6312_image > flux6312_levels[i])} voxels')

    for i in np.arange(len(flux6563_levels)):
        print(f'{i}) H1_6563A band image {invers_pertil_array[i]} percentil : flux {flux6563_levels[i]:.1f} for {np.sum(flux6563_image > flux6563_levels[i])} voxels')

    # Loop throught the regions
    region_dict = {}
    for idx_contour in np.arange(Halpha_idxMinPercentil):

        # Search within that limit
        if idx_contour < SIII_idxMinPercentil:

            if idx_contour == 0:
                maFlux_image = np.ma.masked_where((flux6312_image >= flux6312_levels[idx_contour]) &
                                                  (flux6563_image > Halpha_min_level),
                                                  flux6563_image)
            else:
                maFlux_image = np.ma.masked_where((flux6312_image >= flux6312_levels[idx_contour]) &
                                                  (flux6312_image < flux6312_levels[idx_contour - 1]) &
                                                  (flux6563_image > Halpha_min_level),
                                                  flux6563_image)
        elif idx_contour == 3:
            maFlux_image = np.ma.masked_where((flux6563_image > flux6563_levels[idx_contour]) &
                                              (flux6312_image < flux6312_levels[idx_contour - 1]) &
                                              (flux6563_image > Halpha_min_level),
                                              flux6563_image)

        else:
            maFlux_image = np.ma.masked_where((flux6563_image >= flux6563_levels[idx_contour]) &
                                              (flux6563_image < flux6563_levels[idx_contour-1]) &
                                              (flux6563_image > Halpha_min_level),
                                              flux6563_image)

        region_dict[f'mask_{idx_contour}'] = maFlux_image


    for idx_region, region_items in enumerate(region_dict.items()):
        region_label, region_mask = region_items

        message = f'Mask {region_label} with {np.sum(region_mask.mask)} voxels ' \
                  f'with minimum S3_6312A flux: {np.min(flux6312_image[region_mask.mask]):.1f} and ' \
                  f'with minimum H1_6563A flux: {np.min(flux6563_image[region_mask.mask]):.1f}'
        print(message)



    # # Plot combined mask
    # defaultConf = STANDARD_PLOT.copy()
    # rcParams.update(defaultConf)
    #
    # fig = plt.figure(figsize=(12, 8))
    # ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
    # im = ax.imshow(flux6563_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=flux6563_levels[-2],
    #                                                                     vmin=flux6563_levels[-2],
    #                                                                     base=10))
    #
    # halpha_cmap = cm.gray
    #
    # cmap = cm.get_cmap('viridis', len(region_dict))
    # legend_list = [None] * len(region_dict)
    # alpha_levels = np.linspace(0.1, 0.5, len(region_dict))[::-1]
    #
    # for idx_region, region_items in enumerate(region_dict.items()):
    #
    #     region_label, region_mask = region_items
    #     inv_mask_array = np.ma.masked_array(region_mask.data, ~region_mask.mask)
    #     print(region_label, np.sum(region_mask.mask), flux6312_levels[idx_region])
    #
    #     cm_i = colors.ListedColormap(['black', cmap(idx_region)])
    #     legend_list[idx_region] = patches.Patch(color=cmap(idx_region), label=f'Mask: {region_label}')
    #
    #     ax.imshow(inv_mask_array, cmap=cm_i, vmin=0, vmax=1, alpha=alpha_levels[idx_region])
    #
    # ax.legend(handles=legend_list,  bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    # ax.update({'title': r'{} masks'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
    # ax.set_xlim(90, 220)
    # ax.set_ylim(70, 240)
    # plt.show()
    # plt.savefig(masks_plot)

    # # Store the mask
    # hdul_masks = fits.HDUList()
    # hdul_masks.append(fits.PrimaryHDU())
    # for idx_region, region_items in enumerate(region_dict.items()):
    #     region_label, region_mask = region_items
    #     mask_hdu = fits.ImageHDU(name=region_label, data=region_mask.mask.astype(int), ver=1)
    #     hdul_masks.append(mask_hdu)
    # hdul_masks.writeto(masks_fits, overwrite=True, output_verify='fix')
