import numpy as np
import src.specsiser as sr
import lineid_plot
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, colors, cm, ticker, gridspec
from matplotlib.cbook import get_sample_data
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from lime.tools import label_decomposition

def double_line_blended(lineList):
    new_list = []
    for idx, line in enumerate(lineList):
        if '+' in line:
            new_list.append(line.replace('+', '$\n+\n$'))
        else:
            new_list.append(line)
    return new_list

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])
plots_folder = Path(obsData['file_information']['images_folder'])

STANDARD_PLOT = {'figure.figsize': (14, 18),
                 'axes.titlesize': 14,
                 'axes.labelsize': 14,
                 'legend.fontsize': 12,
                 'xtick.labelsize': 12,
                 'ytick.labelsize': 12}
rcParams.update(STANDARD_PLOT)
format_image = {}
format_image['OffsetImage'] = {'zoom': 0.50}
format_image['AnnotationBbox'] = {'xy': (0.89, 0.65),
                                  'xybox': (0., 0.),
                                  'xycoords': 'axes fraction',
                                  'boxcoords': "offset points", "pad": 0.01}

ak_big = lineid_plot.initial_annotate_kwargs()
ak_big['arrowprops']['relpos'] = (0.5, 0.0)
ak_big['rotation'] = 90

ak_small = lineid_plot.initial_annotate_kwargs()
ak_small['arrowprops'] = dict(width=1e-20, headwidth=1e-20, headlength=1e-20, shrink=1e-20, lw=0.01)
# ak_small['arrowprops']['relpos'] = (0.5, 0.0)
ak_small['rotation'] = 90

pk = lineid_plot.initial_plot_kwargs()
pk['linewidth'] = 0.5

a = dict(width=0.1, headwidth=0.1, headlength=0.1, shrink=None)

labels_dict = {'gp030321': ['O2_3726A_m', 'He1_4026A', 'O3_4363A', 'Ar4_4740A', 'O3_5007A', 'H1_6563A_b', 'He1_7065A'],
               'gp101157': ['Ne3_3869A', 'S2_4069A', 'He1_4471A', 'H1_4861A', 'He1_5876A', 'He1_6678A', 'O2_7319A_m'],
               'gp121903': ['H1_3889A_m', 'H1_4102A', 'He2_4686A', 'O3_4959A', 'O1_6300A_b', 'S2_6716A_b']}

plot_x_low, plot_x_high = 3450, 7650

fig = plt.figure(dpi=300)
gs = fig.add_gridspec(nrows=3, ncols=1)

counter = 0
for i, obj in enumerate(['gp030321', 'gp101157', 'gp121903']):

    z = z_objs[i]
    wmin, wmax = wmin_array[i], wmax_array[i]
    fit_conf = obsData[f'{obj}_line_fitting']

    for ext in ['_BR']:

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objFolder = resultsFolder/f'{obj}'
        objMask = objFolder/f'{obj}{ext}_mask.txt'
        lineLog_file, lineGrid_file = objFolder/f'{obj}{ext}_linesLog.txt', objFolder/f'{obj}{ext}_lineGrid.png'
        pdfTableFile, txtTableFile = objFolder/f'{obj}{ext}_linesTable', objFolder/f'{obj}{ext}_linesTable.txt'
        image_obj = objFolder/f'{obj}_SDSS.png'
        print(f'\n-- Treating {i} :{obj}{ext}.fits')

        gs_obj = gs[i].subgridspec(2, 1, height_ratios=[2.5, 1], hspace=0.0, wspace=0.0)
        ax_big = fig.add_subplot(gs_obj[0, :])
        ax_small = fig.add_subplot(gs_obj[1, :])

        if obj in ['gp030321', 'gp101157', 'gp121903']:

            # Load spectrum
            wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
            flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
            lm = sr.LineMesurer(wave, flux, redshift=z, normFlux=flux_norm, crop_waves=(wmin, wmax))

            # Labels for object
            obj_labels = labels_dict[obj]
            ion_array, wave_array, latex_array = label_decomposition(input_lines=obj_labels, combined_dict=obsData[f'{obj}_line_fitting'])
            latex_array = double_line_blended(latex_array)

            # Big spectrum
            ax_big.step(lm.wave_rest, lm.flux, color='tab:blue', label=obj.replace('gp','GP'), linewidth=1)
            ax_big.xaxis.set_major_locator(plt.NullLocator())
            ax_big.xaxis.set_ticklabels([])
            if 2*i == 2:
                ax_big.set_ylabel(r'Flux $(10^{-14}\,erg\,cm^{-2} s^{-1} \AA^{-1})$')
            ax_big.legend(loc=2)
            ax_big.set_xlim(plot_x_low, plot_x_high)

            lineid_plot.plot_line_ids(lm.wave_rest, lm.flux, wave_array, latex_array, ax=ax_big,
                                      annotate_kwargs=ak_big, plot_kwargs=pk, label1_size=9)

            # Small spectrum
            ax_small.step(lm.wave_rest, lm.flux, color='tab:blue', linewidth=0.5)
            low_limit, up_limit = np.median(lm.flux)/2, np.median(lm.flux)*5
            ax_small.set_ylim(low_limit, up_limit)
            ax_small.set_yscale('log')
            ax_small.xaxis.set_major_locator(plt.NullLocator())
            ax_small.yaxis.set_major_locator(plt.NullLocator())
            ax_small.yaxis.set_ticklabels([], minor=True)

            if 2*i+1 == 5:
                ax_small.xaxis.set_major_locator(ticker.FixedLocator([3500, 4500, 5500, 6500, 7500]))
                ax_small.set_xlabel(r'Wavelength $(\AA)$')
            ax_small.set_xlim(plot_x_low, plot_x_high)
            lineid_plot.plot_line_ids(lm.wave_rest, lm.flux, wave_array, ['']*len(wave_array), ax=ax_small,
                                      annotate_kwargs=ak_small, plot_kwargs=pk)

            # Add the image
            fn = get_sample_data(image_obj, asfileobj=False)
            arr_img = plt.imread(fn, format='png')
            imagebox = OffsetImage(arr_img, **format_image['OffsetImage'])
            imagebox.image.axes = ax_big
            ab = AnnotationBbox(imagebox, **format_image['AnnotationBbox'])
            ax_big.add_artist(ab)

plt.tight_layout()
plt.savefig(plots_folder/'sample_spectra_1200DPI.png')
# plt.show()




# import numpy as np

# import pandas as pd
# import src.specsiser as sr
# from pathlib import Path
# from matplotlib import pyplot as plt, rcParams, colors, cm, ticker, gridspec
# from matplotlib.cbook import get_sample_data
# from matplotlib.offsetbox   import OffsetImage, AnnotationBbox
#
# conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
# obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)
#
# objList = obsData['file_information']['object_list']
# dataFolder = Path(obsData['file_information']['data_folder'])
# resultsFolder = Path(obsData['file_information']['results_folder'])
# fileList = obsData['file_information']['files_list']
#
# z_objs = obsData['sample_data']['z_array']
# wmin_array = obsData['sample_data']['wmin_array']
# wmax_array = obsData['sample_data']['wmax_array']
# flux_norm = obsData['sample_data']['norm_flux']
# noise_region = obsData['sample_data']['noiseRegion_array']
# idx_band = int(obsData['file_information']['band_flux'])
# plots_folder = Path('D:\Dropbox\Astrophysics\Papers\gtc_greenpeas\Greenpeas_osiris_manuscript_v0\images')
#
# STANDARD_PLOT = {'figure.figsize': (7, 20),
#                  'axes.titlesize': 14,
#                  'axes.labelsize': 14,
#                  'legend.fontsize': 12,
#                  'xtick.labelsize': 12,
#                  'ytick.labelsize': 12}
# rcParams.update(STANDARD_PLOT)
# format_image = {}
# format_image['OffsetImage'] = {'zoom': 0.30}
# format_image['AnnotationBbox'] = {'xy': (0.89, 0.65),
#                                   'xybox': (0., 0.),
#                                   'xycoords': 'axes fraction',
#                                   'boxcoords': "offset points", "pad": 0.01}
#
#
# plot_x_low, plot_x_high = 3450, 7650
#
# # fig = plt.figure(constrained_layout=False)
# # gs = gridspec.GridSpec(6, 1, height_ratios=[2.5, 1, 2.5, 1, 2.5, 1])
#
# fig = plt.figure(constrained_layout=True)
# gs = fig.add_gridspec(nrows=3, ncols=1)
#
# counter = 0
# for i, obj in enumerate(['gp030321', 'gp101157', 'gp121903']):
#
#     z = z_objs[i]
#     wmin, wmax = wmin_array[i], wmax_array[i]
#     fit_conf = obsData[f'{obj}_line_fitting']
#
#     for ext in ['_BR']:
#
#         # Declare files location
#         fits_file = dataFolder/f'{obj}{ext}.fits'
#         objFolder = resultsFolder/f'{obj}'
#         objMask = objFolder/f'{obj}{ext}_mask.txt'
#         lineLog_file, lineGrid_file = objFolder/f'{obj}{ext}_linesLog.txt', objFolder/f'{obj}{ext}_lineGrid.png'
#         pdfTableFile, txtTableFile = objFolder/f'{obj}{ext}_linesTable', objFolder/f'{obj}{ext}_linesTable.txt'
#         image_obj = objFolder/f'{obj}_SDSS.png'
#         print(f'\n-- Treating {i} :{obj}{ext}.fits')
#
#         gs_obj = gs[i].subgridspec(2, 1, height_ratios=[2.5, 1], hspace=0.0, wspace=0.0)
#         ax_big = fig.add_subplot(gs_obj[0, :])
#         ax_small = fig.add_subplot(gs_obj[1, :])
#
#         if obj in ['gp030321', 'gp101157', 'gp121903']:
#             print('- Ax big', 2*i)
#             print('- Ax small', 2*i+1)
#             # ax_big = fig.add_subplot(gs[2*i, :])
#             # ax_small = fig.add_subplot(gs[2*i+1, :])
#
#             # Load spectrum
#             wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
#             flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
#             lm = sr.LineMesurer(wave, flux, redshift=z, normFlux=flux_norm, crop_waves=(wmin, wmax))
#
#             # Big spectrum
#             ax_big.step(lm.wave, lm.flux, color='tab:blue', label=obj.replace('gp','GP'), linewidth=1)
#             # ax_big.set_xlim(plot_x_low, plot_x_high)
#             ax_big.xaxis.set_major_locator(plt.NullLocator())
#             ax_big.xaxis.set_ticklabels([])
#             if 2*i == 2:
#                 ax_big.set_ylabel(r'Flux $(10^{-14}\,erg\,cm^{-2} s^{-1} \AA^{-1})$')
#             ax_big.legend(loc=2)
#
#
#             # Small spectrum
#             ax_small.step(lm.wave, lm.flux, color='tab:blue', linewidth=0.5)
#             low_limit, up_limit = np.median(lm.flux)/2, np.median(lm.flux)*5
#             ax_small.set_ylim(low_limit, up_limit)
#             # ax_small.set_xlim(plot_x_low, plot_x_high)
#             ax_small.set_yscale('log')
#             ax_small.xaxis.set_major_locator(plt.NullLocator())
#             ax_small.yaxis.set_major_locator(plt.NullLocator())
#             # ax_small.yaxis.set_minor_locator(plt.NullLocator())
#             ax_small.yaxis.set_ticklabels([], minor=True)
#             # ax_small.minorticks_off()
#             # ax_small.xaxis.set_major_locator(ticker.NullFormatter())
#             # ax_small.xaxis.set_minor_formatter(ticker.NullFormatter())
#
#             if 2*i+1 == 5:
#                 ax_small.xaxis.set_major_locator(ticker.FixedLocator([3500, 4500, 5500, 6500, 7500]))
#                 ax_small.set_xlabel(r'Wavelength $(\AA)$')
#
#             # Add the image
#             fn = get_sample_data(image_obj, asfileobj=False)
#             arr_img = plt.imread(fn, format='png')
#             imagebox = OffsetImage(arr_img, **format_image['OffsetImage'])
#             imagebox.image.axes = ax_big
#             ab = AnnotationBbox(imagebox, **format_image['AnnotationBbox'])
#             ax_big.add_artist(ab)
#
# plt.tight_layout()
# # plt.savefig(plots_folder/f'sample_spectra2.png', pdi=300)
# plt.show()