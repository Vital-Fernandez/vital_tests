import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from src.specsiser.print.plot import STANDARD_PLOT
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, image_array_binning
from astropy.io import fits

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
coordinates_keys_list = obsData['data_location']['wcs_key_list']

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']

# Plot set up
defaultConf = STANDARD_PLOT.copy()
rcParams.update(defaultConf)

# Data location
i = 0
obj = objList[i]
cube_address_i = fitsFolder/fileList[i]
objFolder = resultsFolder/obj

# Load data
wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')
wave = wave / (1 + z_objs[i])
print(f'\n- {obj}: Cube dimensions {cube.shape}')

lineLabel = 'H1_6563A'
lineLimits = lineAreas[lineLabel]

line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
flux_image = line_image.data.data
levelContours = np.nanpercentile(flux_image, pertil_array)

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
ax.update({'title': r'{} galaxy, $H\alpha$ flux'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
im = ax.imshow(flux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=levelContours[1], vmin=levelContours[1], base=10))
cntr1 = ax.contour(flux_image, levels=[levelContours[2]], colors='yellow', alpha=0.1)
h1, l1 = cntr1.legend_elements()
legend = ax.legend([h1[0]], [f'$P_{{{pertil_array[2]}}}$ $F_{{\lambda}}$'], facecolor="white")
plt.show()
