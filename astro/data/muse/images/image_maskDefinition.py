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
flux_H1_6563A = line_image.data.data

lineLabel = 'S3_6312A'
lineLimits = lineAreas[lineLabel]
line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
flux_S3_6312A = line_image.data.data

H1_6563A_Contours = np.nanpercentile(flux_H1_6563A, pertil_array)
S3_6312A_Contours = np.nanpercentile(flux_S3_6312A, pertil_array)

sulfur_bdry = 4
hydrogen_bdry = 2

maFlux_image = np.ma.masked_where((flux_S3_6312A >= S3_6312A_Contours[sulfur_bdry]) &
                                  (flux_S3_6312A < S3_6312A_Contours[sulfur_bdry+1]),
                                  flux_H1_6563A)

fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
ax.update({'title': r'{} galaxy, $H\alpha$ flux'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
im = ax.imshow(maFlux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=H1_6563A_Contours[1], vmin=H1_6563A_Contours[1], base=10))
cntr1 = ax.contour(flux_S3_6312A, levels=S3_6312A_Contours[sulfur_bdry:], cmap='viridis', norm=colors.LogNorm())
cntr2 = ax.contour(flux_H1_6563A, levels=[H1_6563A_Contours[hydrogen_bdry]], colors='red', alpha=0.5)

# for idx, percentile in enumerate(pertil_array[sulfur_bdry:]):
#     label = r'$P_{{{}}}([SIII]6312\AA)$'.format(percentile)
#     cntr1.collections[idx].set_label(label)
# cntr2.collections[0].set_label(r'$P_{{{}}}(H\alpha)$'.format(pertil_array[hydrogen_bdry]))
# ax.legend()

plt.show()
