import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, STANDARD_PLOT, store_frame_to_fits
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
import time

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
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

# Plot set up
defaultConf = STANDARD_PLOT.copy()
rcParams.update(defaultConf)

# Data location
i = 0
obj = objList[i]
cube_address = fitsFolder / fileList[i]
objFolder = resultsFolder / obj
voxelFolder = resultsFolder / obj / 'voxel_data'
db_addresss = objFolder / f'{obj}_database.fits'
mask_address = dataFolder / obj / f'{obj}_mask.txt'

# Declare voxels to analyse
flux6312_image = fits.getdata(db_addresss, f'{ref_flux_line}_flux', ver=1)
flux6312_levels = np.nanpercentile(flux6312_image, pertil_array)

flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
flux6563_levels = np.nanpercentile(flux6563_image, pertil_array)
hdr_plot = fits.getheader(db_addresss, extname='PlotConf')

region_idcs = np.array([0, 1, 2, 3, 4])
for idx in region_idcs:

    mask_label = f'region_{idx}'
    mask_array = fits.getdata(db_addresss, mask_label, ver=1)

    flux_Image = np.ma.masked_array(flux6563_image, mask=mask_array)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))
    ax.update({'title': r'{} galaxy, {}, $H\alpha$ flux'.format(obj, mask_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    im = ax.imshow(flux_Image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=flux6563_levels[-2], vmin=flux6563_levels[-2], base=10))

    plt.show()

