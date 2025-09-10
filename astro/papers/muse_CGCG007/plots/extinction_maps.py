import lime
import numpy as np
import lime as lm
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors, gridspec
from astropy.wcs import WCS
from astropy.io import fits
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits
from mpl_toolkits.axes_grid1 import ImageGrid
import pyneb as pn

# Store emissivity ratios at standard conditions
H1 = pn.RecAtom('H', 1)
temp, den = 10000.0, 100.0
theoEmis_dict = {}

lines_dict = {'H1_6563A': r'$\frac{H\alpha}{H\beta}$',
              'H1_9229A': r'$\frac{H_{Paschen, 9}}{H\beta}$',
              'H1_9015A': r'$\frac{H_{Paschen, 10}}{H\beta}$',
              'H1_8863A': r'$\frac{H_{Paschen, 11}}{H\beta}$'}
max_values = [5.0, 0.05, 0.05, 0.05, 0.05]

for chemLabel, plotLabel in lines_dict.items():
    ion, wave, latexLabel = lm.label_decomposition(chemLabel, scalar_output=True)
    theoRatio = H1.getEmissivity(temp, den, wave=wave) / H1.getEmissivity(temp, den, wave=4861)
    theoEmis_dict[chemLabel] = theoRatio

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

STANDARD_PLOT = {'figure.figsize': (20, 14), 'axes.titlesize': 16,
                 'axes.labelsize': 14, 'legend.fontsize': 10,
                 'xtick.labelsize': 10, 'ytick.labelsize': 10}
rcParams.update(STANDARD_PLOT)

i, obj = 0, 'CGCG007'
lineLabel, lineLimits = 'H1_6563A', (6558.0, 6568.0)


lines_list = list(lines_dict.keys())
ion_array, wave_array, latex_array = lime.label_decomposition(lines_list)

# Data location
cube_address_i = fitsFolder/fileList[i]
objFolder = resultsFolder/obj
db_addresss = objFolder/f'{obj}_database.fits'

# Load data
wave, cube, header = import_muse_fits(cube_address_i)
wave = wave / (1 + z_objs[i])

# Extract cube slice using mpdaf defult pyPopstar.
line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
flux_image = line_image.data.data
levelContours = np.nanpercentile(flux_image, pertil_array)
min_background_percentil = levelContours[2]
norm_bg = colors.SymLogNorm(linthresh=min_background_percentil, vmin=min_background_percentil, base=10)

wcs = WCS(cube.data_header)
fig = plt.figure(figsize=(10, 10))
gs1 = gridspec.GridSpec(ncols=2, nrows=2, width_ratios=[1, 1])
gs1.update(wspace=0.01, hspace=0.02)

parameter_fits = objFolder / 'gauss_flux.fits'
Hbeta_flux = fits.getdata(parameter_fits, 'H1_4861A')

cmap = plt.cm.get_cmap("RdBu_r")
# cmap.set_bad(color='black')

# cmap.set_under("magenta")
# cmap.set_over("yellow")

for i, line in enumerate(lines_list):

    flux_line = fits.getdata(parameter_fits, line)
    ax = fig.add_subplot(gs1[i])#, projection=wcs[0])

    coeff_flux = flux_line / Hbeta_flux
    divnorm = colors.TwoSlopeNorm(vmin=0.0, vcenter=theoEmis_dict[line], vmax=max_values[i])

    ax.imshow(flux_image, cmap=cm.gray, norm=norm_bg, interpolation='none')
    ax.imshow(coeff_flux, cmap='RdBu_r', norm=divnorm, interpolation='none')

    # Text with the line name
    trans = ax.get_xaxis_transform()
    ax.text(0.05, 0.80, lines_dict[line], transform=ax.transAxes, fontsize=35, color='black')

    ax.yaxis.set_major_locator(plt.NullLocator())
    ax.yaxis.set_ticklabels([])
    ax.xaxis.set_major_locator(plt.NullLocator())
    ax.xaxis.set_ticklabels([])
    ax.update({'xlabel': ' '})
    ax.set_ylabel(' ')

    # ax.set_xlim(95, 205)
    # ax.set_ylim(75, 225)
    ax.set_xlim(150, 180)
    ax.set_ylim(150, 180)

# plt.show()
output_file = plotsFolder/f'extinction_maps.png'
print(output_file)
plt.savefig(output_file, bbox_inches='tight')




# import matplotlib.pyplot as plt
#
# from astropy.wcs import WCS
# from astropy.io import fits
# from astropy.utils.data import get_pkg_data_filename
#
# filename = get_pkg_data_filename('galactic_center/gc_msx_e.fits')
#
# hdu = fits.open(filename)[0]
# wcs = WCS(hdu.header)
#
# ax = plt.subplot(projection=wcs[0])
# ax.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
# ax.grid(color='white', ls='solid')
# ax.update({'xlabel': 'Galactic Longitude', 'ylabel': ''})
# ax.yaxis.set_ticklabels([], minor=True)
# ax.yaxis.set_major_locator(plt.NullLocator())
# plt.show()


# # Plot the image:
# wcs = WCS(cube.data_header)
# fig = plt.figure(figsize=(12, 8))#, dpi=600)
# ax = fig.add_subplot(projection=wcs[0])
#
# normalization_background = colors.SymLogNorm(linthresh=min_background_percentil, vmin=min_background_percentil, base=10)
# im = ax.imshow(flux_image, cmap=cm.gray, norm=normalization_background, interpolation='none')
#
# # # Plot contours
# # min_percentil_background = 1
# # contours_levels = levelContours[min_percentil_background:]
# # cntr1 = ax.contour(flux_image, levels=contours_levels, cmap='viridis', norm=colors.LogNorm())
# # for idx, percentile in enumerate(pertil_array[min_percentil_background:]):
# #     label = r'$P_{{{}}}$({})'.format(pertil_array[idx], latexLabel)
# #     cntr1.collections[idx].set_label(label)
# # ax.legend()
#
# ax.set_xlim(95, 205)
# ax.set_ylim(75, 225)
# ax.xaxis.set_major_locator(plt.NullLocator())
# ax.yaxis.set_major_locator(plt.NullLocator())
# ax.yaxis.set_ticklabels([], minor=True)
# ax.update({'xlabel': r'RA', 'ylabel': r'DEC'})
# plt.show()
# # # plt.savefig(plotsFolder/'CGCG007_halpha_image_noAxis.png', bbox_inches='tight')
