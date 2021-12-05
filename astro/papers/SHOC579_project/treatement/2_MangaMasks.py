import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors, patches
from astropy.wcs import WCS
from astropy.io import fits
import lime as lm
from astro.papers.SHOC579_project.SHOC579_methods import open_manga_cubes, DARK_PLOT, line_regions

obs_conf = lm.load_cfg(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\obsConf.ini')

files_data = obs_conf['data_location_windows']
data_folder = Path(files_data['data_folder'])
results_folder = Path(files_data['results_folder'])
file_list = files_data['file_list']
obj_list = files_data['obj_list']

z_objs = obs_conf['sample_data']['z_array']
percentil_array = obs_conf['sample_data']['percentil_array']
percentil_Halpha_array = obs_conf['sample_data']['percentil_Halpha_array']
coordinates_keys_list = obs_conf['sample_data']['wcs_key_list']

dict_errs = {}
dict_nan_values = {}

ref_flux_line = 'O3_4363A'

for i, obj in enumerate(obj_list):

    # Data location
    cube_address_i = data_folder/file_list[i]
    objFolder = results_folder/obj
    db_address = objFolder/f'{obj}_database.fits'

    # Output data:
    hdul_lineslog = fits.HDUList()
    masks_fits = objFolder/f'{obj}_masks.fits'
    masks_plot = objFolder/f'{obj}_masks.png'

    # Load data
    wave, flux, err, hdr = open_manga_cubes(cube_address_i)
    print(f'\n- {obj}: Cube dimensions {flux.shape}')

    # Extract cube slice using mpdaf defult tools.
    idcs_line = np.searchsorted(wave, line_regions[ref_flux_line])
    flux4363_image = flux[idcs_line[0]:idcs_line[1], :, :].sum(axis=0)
    flux4363_levels = np.nanpercentile(flux4363_image, percentil_array)
    flux4363_levels = flux4363_levels[::-1]

    flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    flux6563_levels = np.nanpercentile(flux6563_image, percentil_Halpha_array)

    # fig = plt.figure(figsize=(12, 8))
    # ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y', 1))
    # ax.contour(flux4363_image, levels=flux4363_levels, cmap='viridis', norm=colors.LogNorm())
    # im = ax.imshow(np.log10(flux6563_image), cmap=cm.gray)
    # plt.show()

    # Looping in inverse order
    region_dict = {}
    counter_region = 0
    for idx_contour in np.arange(0, 5):

        if idx_contour == 0:
            maFlux_image = np.ma.masked_where(flux4363_image >= flux4363_levels[idx_contour],
                                              flux4363_image)
        elif idx_contour < 4:
            maFlux_image = np.ma.masked_where((flux4363_image >= flux4363_levels[idx_contour]) &
                                              (flux4363_image < flux4363_levels[idx_contour - 1]),
                                              flux4363_image)
        else:
            maFlux_image = np.ma.masked_where((flux4363_image >=flux4363_levels[8]) &
                                              (flux4363_image < flux4363_levels[idx_contour-1]),
                                              flux4363_image)

        maFlux_image.mask[:20, :] = False
        region_dict[f'region_{idx_contour}'] = maFlux_image

        # fig = plt.figure(figsize=(12, 8))
        # ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y', 1))
        # ax.contour(flux4363_image, levels=flux4363_levels[::-1], cmap='viridis', norm=colors.LogNorm())
        # im = ax.imshow(np.log10(maFlux_image), cmap=cm.gray)
        # plt.show()


    # Plot combined mask
    defaultConf = DARK_PLOT.copy()
    rcParams.update(defaultConf)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y', 1))
    im = ax.imshow(flux6563_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=flux6563_levels[1],
                                                                        vmin=flux6563_levels[1],
                                                                        base=10))

    cmap = cm.get_cmap('viridis', len(region_dict))
    legend_list = [None] * len(region_dict)
    alpha_levels = np.linspace(0.1, 0.5, len(region_dict))[::-1]

    for idx_region, region_items in enumerate(region_dict.items()):

        region_label, region_mask = region_items
        print(f'{region_label}', np.sum(region_mask.mask))
        inv_mask_array = np.ma.masked_array(region_mask.data, ~region_mask.mask)

        cm_i = colors.ListedColormap([cmap(idx_region)])
        legend_list[idx_region] = patches.Patch(color=cmap(idx_region), label=f'Mask: {region_label}')

        ax.imshow(inv_mask_array, cmap=cm_i, vmin=0, vmax=1, alpha=alpha_levels[idx_region])

    # ax.contour(flux4363_image, levels=flux4363_levels[::-1][2:], cmap='viridis_r', norm=colors.LogNorm())
    ax.legend(handles=legend_list,  bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    ax.update({'title': r'{} masks'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
    plt.show()
    plt.savefig(masks_plot)

    # Store the mask
    hdul_masks = fits.HDUList()
    hdul_masks.append(fits.PrimaryHDU())
    for idx_region, region_items in enumerate(region_dict.items()):
        region_label, region_mask = region_items
        mask_hdu = fits.ImageHDU(name=region_label, data=region_mask.mask.astype(int), ver=1)
        hdul_masks.append(mask_hdu)
    hdul_masks.writeto(masks_fits, overwrite=True, output_verify='fix')
