import numpy as np
import lime
from lime.plots import STANDARD_PLOT
from src.specsiser.plots import latex_labels
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']

fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

voxel_grid_size = obsData['sample_data']['grid_shape_array']

HII_chim_params = ['OH', 'NO', 'logU']

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    maskFits_address, mask_list = objFolder/f'{obj}_masks.fits', ['MASK_0', 'MASK_1', 'MASK_2']
    db_address = objFolder/f'{obj}_database.fits'
    outputDb = objFolder/f'{obj}_chemical.fits'
    chemFolder = objFolder/'chemistry'
    HII_chim_results = objFolder / f'{obj}_HIIchimistry.fits'

    flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')

    # Data from the direct method
    hdr = fits.getheader(f'{chemFolder}/O2.fits', extname='O2')
    O2_image = fits.getdata(f'{chemFolder}/O2.fits', extname='O2')
    O3_image = fits.getdata(f'{chemFolder}/O3.fits', extname='O3')
    nSII_image = fits.getdata(f'{chemFolder}/n_e.fits', extname='n_e')
    O_direct = O2_image + O3_image


    # Data from HII-CHIMS-TRY the direct method
    # lime.save_param_maps(HII_chim_results, HII_chim_params, lines_list=None, output_folder=chemFolder,
    #                      spatial_mask_file=maskFits_address, ext_mask=mask_list,
    #                      ext_log='_HIICHIMISTRY', row_suffix='_LINE', output_files_prefix='HII_CHIMISTRY_MAP_')

    O_HIICHIMISTRY = fits.getdata(f'{chemFolder}/HII_CHIMISTRY_MAP_OH.fits', ext_name='OH')

    # Plot HII chemistry
    defaultConf = STANDARD_PLOT.copy()
    rcParams.update(defaultConf)

    halpha_cmap = cm.gray
    halpha_cmap.set_under('black')

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y'))

    bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)

    im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
    im2 = ax.imshow(O_HIICHIMISTRY)
    cbar = fig.colorbar(im2, ax=ax)
    param_label = r'$\frac{O}{H}$ direct method'
    ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    ax.set_xlim(120, 210)
    ax.set_ylim(110, 220)
    # plt.savefig(objFolder/f'{obj}_parameter_map_{ext_label}')
    plt.show()

    # # Plot Direct method
    # defaultConf = STANDARD_PLOT.copy()
    # rcParams.update(defaultConf)
    #
    # halpha_cmap = cm.gray
    # halpha_cmap.set_under('black')
    #
    # fig = plt.figure(figsize=(10, 10))
    # ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y'))
    #
    # bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)
    #
    # im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
    # im2 = ax.imshow(O_HIICHIMISTRY)
    # cbar = fig.colorbar(im2, ax=ax)
    # param_label = r'\frac{O}{H} direct method$'
    # ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    # ax.set_xlim(120, 210)
    # ax.set_ylim(110, 220)
    # # plt.savefig(objFolder/f'{obj}_parameter_map_{ext_label}')
    # plt.show()



    # # a['OH'] = 7.709410496198816
    #
