import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import background_color, DARK_PLOT, label_Conver, latex_Conver
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
import pyneb as pn

sulfur_lines = {'S2_6716A': r'$[SII]6716\AA', 'S2_6731A': r'$[SII]6716\AA'}

dinamicLines = {'H1_6563A': r'$H\alpha_{Narrow}$',
                'H1_4861A': r'$H\beta_{Narrow}$',
                'O3_5007A': r'$[OIII]5007\AA_{Narrow}$',
                'O3_5007A': r'$[OIII]5007\AA_{Narrow}$',
                'S3_9069A': r'$[SIII]9069\AA_{Narrow}$',
                'He1_5876A': r'$HeI\,5876\AA_{Narrow}$'}


conf_file = Path('/home/vital/PycharmProjects/vital_tests/astro/data/muse/muse_greenpeas.ini')
seminar_folder = Path('/home/vital/Dropbox/Astrophysics/Seminars/UniVapo 2021/')

# Declare data and files location
obsData = sr.loadConfData(conf_file, group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
voxel_grid_size = obsData['sample_data']['grid_shape_array']
coordinates_keys_list = obsData['data_location']['wcs_key_list']

# Store emissivity ratios at standard conditions
H1 = pn.RecAtom('H', 1)
temp, den = 10000.0, 100.0

theoEmis_dict = {}
for chemLabel, plotLabel in label_Conver.items():
    ion, wave, latexLabel = sr.label_decomposition(chemLabel, scalar_output=True)
    dict_label = f'{plotLabel}/Hbeta'
    theoRatio = H1.getEmissivity(temp, den, wave=wave) / H1.getEmissivity(temp, den, wave=4861)
    theoEmis_dict[dict_label] = theoRatio

i = 0
obj = objList[i]

# Data location
objFolder = resultsFolder/obj
voxelFolder = resultsFolder/obj/'voxel_data'
db_address = objFolder / f'{obj}_database.fits'
# fitsLog_address = objFolder / f'{obj}_linesLog.fits'
chemMaps_fits_address = objFolder / f'{obj}_ChemParamMaps.fits'
theoLineMaps_fits_address = objFolder / f'{obj}_fitTheoFluxesMaps.fits'
linesMaps_fits_address = objFolder/f'{obj}_lineParamMaps.fits'

# ----------------------------------------- Line flux images

hdr_plot = fits.getheader(linesMaps_fits_address, extname='PlotConf')
flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
halpha_min_level = fits.getval(db_address, keyword=f'P9000', extname=f'H1_6563A_flux')
halpha_thresd_level = fits.getval(db_address, keyword=f'P8000', extname=f'H1_6563A_flux')
# linthresh = flux6563_levels[-2], vmin = flux6563_levels[-3],

defaultConf = DARK_PLOT.copy()
rcParams.update(defaultConf)

halpha_cmap = cm.gray
halpha_cmap.set_under(background_color)

# # Recombination lines
# for chemLabel, plotLabel in label_Conver.items():
#     fig = plt.figure(figsize=(10, 10))
#     ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))
#
#     dict_label = f'{plotLabel}/Hbeta'
#     flux_image = fits.getdata(linesMaps_fits_address, dict_label, ver=1)
#
#     divnorm = colors.TwoSlopeNorm(vmin=np.nanmin(flux_image),
#                                   vcenter=theoEmis_dict[dict_label],
#                                   vmax=np.nanmax(flux_image))
#
#     im = ax.imshow(flux6563_image, cmap=halpha_cmap,
#                    norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))
#     im2 = ax.imshow(flux_image, cmap='RdBu', norm=divnorm)
#
#     cbar = fig.colorbar(im2, ax=ax)
#     cbar.set_label('Line ratio (white theoretical value)', rotation=270, labelpad=50, fontsize=15)
#
#     ax.set_xlim(140, 210)
#     ax.set_ylim(110, 240)
#
#     ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_Conver[chemLabel], latex_Conver['H1_4861A'])
#     ax.update({'title': r'Galaxy {}: {}'.format(obj, ratio_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
#     plt.savefig(seminar_folder / f'{chemLabel}_recomb', dpi=300, bbox_inches='tight')

    # plt.tight_layout()
    # plt.show()

dinamicLines = {'H1_6563A': r'$H\alpha_{Narrow}$',
                'H1_4861A': r'$H\beta_{Narrow}$',
                'O3_5007A': r'$[OIII]5007\AA_{Narrow}$'}
# dinamicLines.update(sulfur_lines)
('v_r', 'sigma_vel')
# Recombination lines
for dinLabel, latex_label in dinamicLines.items():
for param in ['v_r']:

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot()

        dict_label = f'{param}_{dinLabel}'
        param_image = fits.getdata(linesMaps_fits_address, dict_label, ver=1)

        im = ax.imshow(flux6563_image, cmap=halpha_cmap,
                       norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))

        ax.set_xlim(100, 210)
        ax.set_ylim(90, 240)

        if param == 'v_r':
            # param_min, param_max = np.nanmin(param_image), np.nanmax(param_image)
            # divnorm = colors.TwoSlopeNorm(vmin=param_min,
            #                               vcenter=0.0,
            #                               vmax=param_max)
            # im2 = ax.imshow(param_image, cmap='RdBu', norm=divnorm)
            im2 = ax.imshow(param_image, vmin=-5, vmax=80)

        else:
            im2 = ax.imshow(param_image, vmin=30, vmax=80)

        cbar = fig.colorbar(im2)
        label_bar = latex_Conver[param]
        cbar.set_label(label_bar, rotation=270, labelpad=50, fontsize=20)

        ax.update({'title': r'Galaxy {}: {}'.format(obj, latex_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
        # plt.savefig(seminar_folder / f'{param}_{dinLabel}_image', dpi=300, bbox_inches='tight')

        plt.tight_layout()
        # plt.savefig(resultsFolder/obj/f'map_{obj}_{dinLabel}_{param}.png', bbox_inches='tight')
        plt.show()


# # ----------------------------------------- Parameter maps
#
# hdr_plot = fits.getheader(chemMaps_fits_address, extname='PlotConf')
# flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
# halpha_min_level = fits.getval(db_address, keyword=f'P9000', extname=f'H1_6563A_flux')
# halpha_thresd_level = fits.getval(db_address, keyword=f'P8000', extname=f'H1_6563A_flux')
#
# defaultConf = DARK_PLOT.copy()
# rcParams.update(defaultConf)
#
# halpha_cmap = cm.gray
# halpha_cmap.set_under(background_color)
#
# # param lines
# for i, param in enumerate(['T_low', 'n_e', 'cHbeta', 'N_O', 'S2', 'S3', 'N2', 'O2', 'O3', 'He1']):
#
#     fig = plt.figure(figsize=(10, 10))
#     ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))
#
#     if param != 'N_O':
#         param_image = fits.getdata(chemMaps_fits_address, param, ver=1)
#     else:
#         N2 = np.power(10, fits.getdata(chemMaps_fits_address, 'N2', ver=1) - 12)
#         O2 = np.power(10, fits.getdata(chemMaps_fits_address, 'O2', ver=1) - 12)
#         param_image = np.log10(N2 / O2)
#         latex_labels['N_O'] = r'log(N/O)'
#
#     ax.set_xlim(140, 210)
#     ax.set_ylim(110, 240)
#
#     im = ax.imshow(flux6563_image, cmap=halpha_cmap,
#                    norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))
#     im2 = ax.imshow(param_image)
#
#     cbar = fig.colorbar(im2, ax=ax)
#
#     param_label = latex_labels[param]
#     ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
#     plt.savefig(seminar_folder/f'{param}_image', dpi=300, bbox_inches='tight')
#     # plt.savefig(plot_address, dpi=300, bbox_inches='tight')
#     # plt.tight_layout()
#
#     # plt.show()

