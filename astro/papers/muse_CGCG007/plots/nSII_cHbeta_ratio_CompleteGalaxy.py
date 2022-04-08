import numpy as np
import src.specsiser as sr
import time
import lime

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import chemical_lines_indexing, import_muse_fits
from astropy.io import fits
from astropy.table import Table
from lime.plots import STANDARD_PLOT
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from src.specsiser.plots import latex_labels
import pyneb as pn


# a = np.array([[1, 2, 3],
#               [4, np.nan, 6],
#               [7, 8, 9]])
# print(a)
#
# idcs_nan = ~np.isnan(a)
# vector = a[idcs_nan].flatten()
#
# print(vector)
#
# a[idcs_nan] = vector * 2
#
# print(a)

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

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

S2 = pn.Atom('S', 2)

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = f'{objFolder}/chemistry'
    obsLog_addresss = objFolder / f'{obj}_linesLog.fits'
    absLog_addresss = objFolder / f'{obj}_linesLog_abs.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    db_address = objFolder / f'{obj}_database.fits'
    ne_SII_map = objFolder/f'{obj}_neSII_gauss.fits'
    extinction_map = objFolder/f'{obj}_HI_extinction.fits'


    neSII_map = fits.getdata(ne_SII_map, extname='ne_SII')
    cHbeta_map = fits.getdata(extinction_map, extname='cHbeta')
    hdr = fits.getheader(ne_SII_map, extname='ne_SII')

    idcs_ne_cHbeta = ~(np.isnan(neSII_map) | np.isnan(cHbeta_map))
    norm_density = neSII_map[idcs_ne_cHbeta]/cHbeta_map[idcs_ne_cHbeta]
    norm_density_map = np.full(neSII_map.shape, np.nan)
    norm_density_map[idcs_ne_cHbeta] = norm_density

    # ----------------------------------------- Generate the image plots ----------------------------------------
    flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')

    defaultConf = STANDARD_PLOT.copy()
    rcParams.update(defaultConf)

    halpha_cmap = cm.gray
    halpha_cmap.set_under('black')

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y'))

    bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)

    im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
    im2 = ax.imshow(norm_density_map, norm=colors.LogNorm())
    cbar = fig.colorbar(im2, ax=ax)
    param_label = latex_labels['n_e']
    ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    ax.set_xlim(95, 210)
    ax.set_ylim(80, 222)
    # plt.savefig(objFolder/f'{obj}_parameter_map_{ext_label}')
    plt.show()


    # Histogram plot
    defaultConf = STANDARD_PLOT.copy()
    rcParams.update(defaultConf)

    fig = plt.figure(figsize=(10, 10))

    ax = fig.add_subplot()

    ax.hist(norm_density_map.flatten(), bins=50, log=True)

    ratio_label = r'$\frac{n_{[SII]}}{c(H\beta)}$'
    voxel_count = np.sum(~np.isnan(norm_density_map))

    ax.update({'xlabel': ratio_label, 'ylabel': r'Voxel count',
               'title': f'{ratio_label} = ${np.nanmedian(norm_density_map):.0f}\pm{np.nanstd(norm_density_map):.0f}$,  direct method ({voxel_count} voxels)'})

    ax.legend()
    plt.show()