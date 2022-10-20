import numpy as np
import pandas as pd
import time
import lime
from pathlib import Path
import lineid_plot
from lime.plots import latex_science_float


from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits
from matplotlib import pyplot as plt, rcParams, ticker
from matplotlib.cbook import get_sample_data
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

ak_big = lineid_plot.initial_annotate_kwargs()
# ak_big['arrowprops']['relpos'] = (0.5, 0.0)
ak_big['rotation'] = 90

ak_small = lineid_plot.initial_annotate_kwargs()
ak_small['arrowprops'] = dict(width=1e-20, headwidth=1e-20, headlength=1e-20, shrink=1e-20, lw=0.01)
# ak_small['arrowprops']['relpos'] = (0.5, 0.0)
ak_small['rotation'] = 90

pk = lineid_plot.initial_plot_kwargs()
pk['linewidth'] = 0.25

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
plotsFolder = Path(obsData['data_location']['plots_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = 1

i = 0
obj = 'CGCG007'
voxel_cords = (167, 167)
lineLabels = np.array(['Ar4_4740A', 'H1_4861A', 'O3_4959A', 'He1_4922A', 'Fe3_4987A', 'O3_5007A', 'N1_5200A', 'Cl3_5518A',
                       'Cl3_5538A', 'He1_5876A', 'O1_6300A', 'S3_6312A', 'O1_6364A', 'H1_6563A', 'N2_6548A', 'N2_6584A',
                       'He1_6678A', 'S2_6716A', 'S2_6731A',
                       'He1_7065A', 'Ar3_7136A', 'He1_7281A', 'O2_7320A', 'O2_7330A', 'S3_9069A', 'H1_8204A', 'H1_9229A',
                       'H1_9015A', 'H1_8863A', 'H1_8750A'])

uniq, count = np.unique(lineLabels, return_counts=True)
np.any(count > 1)

lineIons, lineWaves, lineLatex = lime.label_decomposition(lineLabels)

labels_dict = {'H1_8204A': 'Paschen\nJump', 'H1_6563A': r'$H\alpha$', 'H1_4861A': r'$H\beta$'}
for old_label, new_label in labels_dict.items():
    idx = np.where(lineLabels == old_label)
    lineLatex[idx] = new_label


# Data location
cube_address = fitsFolder/fileList[i]
objFolder = resultsFolder/obj
voxelFolder = resultsFolder/obj/'voxel_data'
db_addresss = objFolder/f'{obj}_database.fits'
maskFits_address = objFolder/f'{obj}_masks.fits'

idx_j, idx_i = voxel_cords
wave, cube, header = import_muse_fits(cube_address)
flux_voxel = cube[:, idx_j, idx_i].data.data / 1000
flux_err = np.sqrt(cube[:, idx_j, idx_i].var.data) / 1000

voxel = lime.Spectrum(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i])
# voxel.plot_spectrum(spec_label=f'{obj} voxel {idx_j}-{idx_i}', log_scale=True)

STANDARD_PLOT = {'figure.figsize': (16, 5),
                 'axes.titlesize': 10,
                 'axes.labelsize': 6,
                 'legend.fontsize': 12,
                 'xtick.labelsize': 6,
                 'ytick.labelsize': 6}
rcParams.update(STANDARD_PLOT)

fig = plt.figure(dpi=600)
gs = fig.add_gridspec(nrows=1, ncols=1)
gs_obj = gs[i].subgridspec(2, 1, height_ratios=[5.5, 1], hspace=0.0, wspace=0.0)
ax_big = fig.add_subplot(gs_obj[0, :])
ax_small = fig.add_subplot(gs_obj[1, :])

# Big spectrum
ax_big.step(voxel.wave_rest, voxel.flux, color='tab:blue', linewidth=0.5)
ax_big.xaxis.set_major_locator(plt.NullLocator())
ax_big.xaxis.set_ticklabels([])
if 2 * i == 2:
    ax_big.set_ylabel(r'Flux $(10^{-14}\,erg\,cm^{-2} s^{-1} \AA^{-1})$')
plot_x_low, plot_x_high = 4700, 9350
ax_big.set_xlim(plot_x_low, plot_x_high)

lineid_plot.plot_line_ids(voxel.wave_rest, voxel.flux, lineWaves, lineLatex, ax=ax_big,
                          annotate_kwargs=ak_big, plot_kwargs=pk, label1_size=6)

plot_labels = {'xlabel': r'Wavelength $(\AA)$',
               'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})$'}
plot_labels['ylabel'] = plot_labels['ylabel'] + r' $\,/\,{}$'.format(latex_science_float(1e-20*1000))
ax_big.update(plot_labels)

# Small spectrum
ax_small.step(voxel.wave_rest, voxel.flux, color='tab:blue', linewidth=0.25)
low_limit, up_limit = np.median(voxel.flux) / 4, np.median(voxel.flux) * 7
ax_small.set_ylim(low_limit, up_limit)
ax_small.set_yscale('log')
# ax_small.xaxis.set_major_locator(plt.NullLocator())
ax_small.yaxis.set_major_locator(plt.NullLocator())
ax_small.yaxis.set_ticklabels([], minor=True)
ax_small.set_xlabel(r'Wavelength $(\AA)$')
ax_small.set_xlim(plot_x_low, plot_x_high)
lineid_plot.plot_line_ids(voxel.wave_rest, voxel.flux, lineWaves, [''] * len(lineWaves), ax=ax_small,
                          annotate_kwargs=ak_small, plot_kwargs=pk)
# ax_small.xaxis.set_major_locator(ticker.FixedLocator([4500, 5500, 6500, 7500, 8500, 9000]))

# Image format
format_image = {}
format_image['OffsetImage'] = {'zoom': 0.040}
format_image['AnnotationBbox'] = {'xy': (0.66, 0.54),
                                  'xybox': (0., 0.),
                                  'xycoords': 'axes fraction',
                                  'boxcoords': "offset points", "pad": 0.01}
# bboxprops=dict(edgecolor='red')


image_obj = plotsFolder/'CGCG007_halpha_image_noAxis.png'
fn = get_sample_data(image_obj, asfileobj=False)
arr_img = plt.imread(fn, format='png')
imagebox = OffsetImage(arr_img, **format_image['OffsetImage'])
imagebox.image.axes = ax_big
ab = AnnotationBbox(imagebox, bboxprops=dict(edgecolor='black', lw=0.25), **format_image['AnnotationBbox'])
ax_big.add_artist(ab)

# plt.tight_layout()
# plt.show()

plot_image_file = plotsFolder/'voxel_spectrum_and_muse_halpha.png'
plt.savefig(plot_image_file, bbox_inches='tight')
