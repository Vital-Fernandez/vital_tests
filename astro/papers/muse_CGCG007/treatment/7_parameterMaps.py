import numpy as np
import lime
from lime.plots import STANDARD_PLOT
from src.specsiser.plots import latex_labels
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import save_log_maps


# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
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
    param_list = np.array(['cHbeta', 'Ar4', 'Ar3', 'O2', 'O3', 'N2', 'He1', 'n_e', 'T_low', 'T_high', 'S2', 'S3'])

    # Data for the astronomical coordinates
    hdr = fits.getheader(maskFits_address, extname='MASK_0')

    # Generate the map files
    save_log_maps(outputDb, param_list, chemFolder, maskFits_address, mask_list, ext_log='_CHEMISTRY_OUTPUTS',
                  page_hdr=hdr)


    # ----------------------------------------- Generate the image plots ----------------------------------------
    flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')

    for param in param_list:

        with fits.open(f'{chemFolder}/{param}.fits') as hdu_list:

            image_data, image_header = hdu_list[param].data, hdu_list[param].header

            print(f'{chemFolder}/{param}.fits', param, np.all(np.isnan(image_data)))

            defaultConf = STANDARD_PLOT.copy()
            rcParams.update(defaultConf)

            halpha_cmap = cm.gray
            halpha_cmap.set_under('black')

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(projection=WCS(image_header), slices=('x', 'y'))

            bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)

            im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
            im2 = ax.imshow(image_data)
            cbar = fig.colorbar(im2, ax=ax)
            param_label = latex_labels[param]
            ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
            ax.set_xlim(120, 210)
            ax.set_ylim(110, 220)
            # plt.savefig(objFolder/f'{obj}_parameter_map_{ext_label}')
            plt.show()
