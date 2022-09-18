import lime
import numpy as np
import lime as lm
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors, gridspec
from astropy.wcs import WCS
from astropy.io import fits
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits, latex_labels
from mpl_toolkits.axes_grid1 import ImageGrid, make_axes_locatable

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
                 'xtick.labelsize': 10, 'ytick.labelsize': 10,
                 'font.family':'Times New Roman'}
rcParams.update(STANDARD_PLOT)

i, obj = 0, 'CGCG007'
lineLabel, lineLimits = 'H1_6563A', (6558.0, 6568.0)

# Data location
cube_address_i = fitsFolder/fileList[i]
objFolder = resultsFolder/obj
db_addresss = objFolder/f'{obj}_database.fits'

# Load data
wave, cube, header = import_muse_fits(cube_address_i)
wave = wave / (1 + z_objs[i])

# Extract cube slice using mpdaf defult tools.
line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
flux_image = line_image.data.data
levelContours = np.nanpercentile(flux_image, pertil_array)
min_background_percentil = levelContours[2]
norm_bg = colors.SymLogNorm(linthresh=min_background_percentil, vmin=min_background_percentil, base=10)

wcs = WCS(cube.data_header)
fig = plt.figure(figsize=(7, 8), dpi=600)
gs1 = gridspec.GridSpec(ncols=3, nrows=3)
gs1.update(wspace=0, hspace=0)

params_list = ['OH', 'NO', 'logU']
direct_method_fits = '/home/vital/Dropbox/Astrophysics/Papers/muse_CGCG007/treatment/CGCG007/CGCG007_chemical.fits'

methodology = 'neural_fitting'
conf_list = ['direct_method', 'minOneErr', 'HIICHImistry-PopStar']
chemistryFolder = objFolder/'chemistry'
files_conv = {'OH': 'logOH', 'NO': 'logNO', 'logU': 'logU', 'eta': 'eta', 'S2_S3': 'S2_S3'}
method_conf = {'direct_method': 'neural_fitting_direct_method',
               'minOneErr': 'GridSampling_minOneErr',
               'HIICHImistry-PopStar': 'HII-CHI-mistry_HIICHImistry-PopStar'}
ext_conf = {'OH': 'logOH', 'NO': 'logNO', 'logU': 'logU', 'eta': 'eta', 'S2_S3': 'S2_S3'}
title_conv = {'direct_method':          'Direct Method',
              'minOneErr':              'Neural sampler,\nOne order uncertainty',
              'HIICHImistry-PopStar':   r'HII-CHI-mistry'}
idx_ax = 0
for param_in in params_list:

    for i, conf in enumerate(conf_list):
        if conf == 'direct_method' and param_in == 'logU':
            param = 'S2_S3'
        else:
            param = param_in
        print(conf, param)


        methodology = method_conf[conf]

        # Color normalization based on measurement distribution
        param_distr = np.array(obsData[methodology][f'{param}_global'])
        median = param_distr[0]
        Per84th, Per16th = param_distr[1], param_distr[2]
        vmin, vmax = median - Per16th * 3, median + Per84th * 3
        divnorm = colors.TwoSlopeNorm(vmin=vmin, vcenter=median, vmax=vmax)

        fits_file = chemistryFolder/f'{conf}_{files_conv[param]}.fits'
        with fits.open(fits_file) as hdul:
            if conf in ['minOneErr']:
                param_ext = ext_conf[param]
            else:
                param_ext = param
            param_image, param_hdr = hdul[param_ext].data, hdul[param_ext].header

        ax = fig.add_subplot(gs1[idx_ax])#, projection=wcs[0])
        im1 = ax.imshow(flux_image, cmap=cm.gray, norm=norm_bg, interpolation='none')
        im2 = ax.imshow(param_image, interpolation='none', vmin=vmin, vmax=vmax)

        ticks_values = np.array([median - Per16th*2, median - Per16th, median, median + Per84th, median + Per84th*2])
        cbar = fig.colorbar(im2, ax=ax, shrink=0.9995, pad=0, orientation='vertical', ticks=ticks_values)
        cbar.ax.yaxis.set_ticks_position('left')
        cbar.ax.tick_params(axis="y", direction="in", colors='white', labelsize=7)

        if idx_ax < 3:
            ax.set_title(title_conv[conf])

        # Text with the line name
        trans = ax.get_xaxis_transform()
        ax.text(0.05, 0.80, latex_labels[param], transform=ax.transAxes, fontsize=16, color='white')

        ax.yaxis.set_major_locator(plt.NullLocator())
        ax.yaxis.set_ticklabels([])
        ax.xaxis.set_major_locator(plt.NullLocator())
        ax.xaxis.set_ticklabels([])
        ax.update({'xlabel': ' '})
        ax.set_ylabel(' ')

        ax.set_xlim(130, 210)
        ax.set_ylim(110, 220)
        idx_ax += 1



# plt.show()
output_file = plotsFolder/f'methodology_maps.png'
print(output_file)
plt.savefig(output_file, bbox_inches='tight')



