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


for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    maskFits_address, mask_list = objFolder/f'{obj}_masks.fits', ['MASK_0', 'MASK_1', 'MASK_2']
    db_address = objFolder/f'{obj}_database.fits'
    outputDb = objFolder/f'{obj}_chemical.fits'
    chemFolder = objFolder/'chemistry'

    flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')

    cHbeta_image = fits.getdata(f'{chemFolder}/cHbeta.fits', extname='cHbeta')
    nSII_image = fits.getdata(f'{chemFolder}/n_e.fits', extname='n_e')
    hdr = fits.getheader(f'{chemFolder}/n_e.fits', extname='n_e')

    param_ratio = nSII_image/cHbeta_image

    defaultConf = STANDARD_PLOT.copy()
    rcParams.update(defaultConf)

    halpha_cmap = cm.gray
    halpha_cmap.set_under('black')

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y'))

    bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)

    im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
    im2 = ax.imshow(param_ratio, norm=colors.LogNorm())
    cbar = fig.colorbar(im2, ax=ax)
    param_label = r'$\frac{n_{[SII]}}{c(H\beta)}$'
    ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    ax.set_xlim(120, 210)
    ax.set_ylim(110, 220)
    # plt.savefig(objFolder/f'{obj}_parameter_map_{ext_label}')
    plt.show()

    # Histogram plot
    defaultConf = STANDARD_PLOT.copy()
    rcParams.update(defaultConf)

    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot()

    ax.hist(param_ratio.flatten(), bins=50, log=True)

    ratio_label = r'$\frac{n_{[SII]}}{c(H\beta)}$'
    voxel_count = np.sum(~np.isnan(param_ratio))

    ax.update({'xlabel': ratio_label, 'ylabel': r'Voxel count',
               'title': f'{ratio_label} = ${np.nanmedian(param_ratio):.0f}\pm{np.nanstd(param_ratio):.0f}$,  direct method ({voxel_count} voxels)'})

    ax.legend()
    plt.show()
