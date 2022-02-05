import numpy as np
import lime as lm
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits

# Declare data and files location
obsData = lm.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
plotsFolder = Path(obsData['data_location']['plots_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']

STANDARD_PLOT = {'figure.figsize': (20, 14), 'axes.titlesize': 16, 'axes.labelsize': 14, 'legend.fontsize': 10,
                 'xtick.labelsize': 10, 'ytick.labelsize': 10}

# Plot set up
# labelsDict = {'xlabel': r'RA',
#               'ylabel': r'DEC'}
rcParams.update(STANDARD_PLOT)

i, obj = 0, 'CGCG007'

# Data location
cube_address_i = fitsFolder/fileList[i]
objFolder = resultsFolder/obj
db_addresss = objFolder/f'{obj}_database.fits'

# Load data
wave, cube, header = import_muse_fits(cube_address_i)
wave = wave / (1 + z_objs[i])
print(f'\n- {obj}: Cube dimensions {cube.shape}')
lineLabel, lineLimits = 'H1_6563A', (6558.0, 6568.0)


plot_image_file = objFolder/f'{obj}_{lineLabel}_contours.png'
ion, wavelength, latexLabel = lm.label_decomposition(lineLabel, scalar_output=True)

# Extract cube slice using mpdaf defult tools.
line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
flux_image = line_image.data.data
levelContours = np.nanpercentile(flux_image, pertil_array)


# Plot the image:
fig = plt.figure(figsize=(10, 10), dpi=600)
# ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
ax = fig.add_subplot()

min_background_percentil = levelContours[2]
normalization_background = colors.SymLogNorm(linthresh=min_background_percentil, vmin=min_background_percentil, base=10)
im = ax.imshow(flux_image, cmap=cm.gray, norm=normalization_background, interpolation='none')

# # Plot contours
# min_percentil_background = 1
# contours_levels = levelContours[min_percentil_background:]
# cntr1 = ax.contour(flux_image, levels=contours_levels, cmap='viridis', norm=colors.LogNorm())
# for idx, percentile in enumerate(pertil_array[min_percentil_background:]):
#     label = r'$P_{{{}}}$({})'.format(pertil_array[idx], latexLabel)
#     cntr1.collections[idx].set_label(label)
# ax.legend()

ax.set_xlim(95, 210)
ax.set_ylim(80, 222)
# ax.update({'title': r'CGCG007-025, {} image'.format(r'$H\alpha$'), 'xlabel': r'RA', 'ylabel': r'DEC'})
ax.xaxis.set_major_locator(plt.NullLocator())
ax.yaxis.set_major_locator(plt.NullLocator())
ax.yaxis.set_ticklabels([], minor=True)

plt.savefig(plotsFolder/'CGCG007_halpha_image_noAxis.png', bbox_inches='tight')
# plt.show()


