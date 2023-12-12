import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import STANDARD_AXES
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from delete.data_printing import DARK_PLOT, background_color, foreground_color

conf_file = Path(r'D:\Pycharm Projects\vital_tests\astro\data\muse\muse_greenpeas.ini')
# conf_file = Path('/home/vital/PycharmProjects/vital_tests/astro/data/muse/muse_greenpeas.ini')

obsData = sr.loadConfData(conf_file, group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

folder = Path('/home/vital/Dropbox/Astrophysics/Seminars/UniVapo 2021/')



# Data location
i = 0
obj = objList[i]
cube_address = fitsFolder / fileList[i]
objFolder = resultsFolder / obj
voxelFolder = resultsFolder / obj / 'voxel_data'
db_addresss = objFolder / f'{obj}_database.fits'
mask_address = dataFolder / obj / f'{obj}_mask.txt'

# Load data
wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')
wave_rest = wave / (1 + z_objs[i])

idx_j, idx_i = 170, 170

flux_voxel = cube[:, idx_j, idx_i].data.data
flux_err = cube[:, idx_j, idx_i].var.data

# Plot set up
defaultConf = DARK_PLOT.copy()

defaultConf['axes.titlesize']= 14
defaultConf['axes.labelsize']= 18
defaultConf['legend.fontsize']= 12
defaultConf['xtick.labelsize']= 14
defaultConf['ytick.labelsize']= 14

rcParams.update(defaultConf)

fig = plt.figure(figsize=(16, 9))
ax = fig.add_subplot()
ax.step(wave_rest, flux_voxel, color=foreground_color)
ax.set_yscale('log')
ax.update({'xlabel': r'Wavelength $(\AA)$',
           'ylabel': r'Flux $(10^{20}\,erg\,cm^{-2} s^{-1} \AA^{-1})$'})

ax.spines['left'].set_position(('outward', 10))
ax.spines['bottom'].set_position(('outward', 10))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()
plt.show()
