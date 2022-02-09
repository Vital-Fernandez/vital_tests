import lime
import pyneb as pn
import numpy as np
import lime as lm

from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS

from astro.papers.muse_CGCG007.muse_CGCG007_methods import label_Conver, target_lines, param_images, latex_Conver, dinamicLines

# Declare data and files location
obsData = lm.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
voxel_grid_size = obsData['sample_data']['grid_shape_array']
coordinates_keys_list = obsData['data_location']['wcs_key_list']

# ------ Extinction coefficients

# # Store emissivity ratios at standard conditions
# H1 = pn.RecAtom('H', 1)
# temp, den = 10000.0, 100.0
# theoEmis_dict = {}
# for chemLabel, plotLabel in label_Conver.items():
#     ion, wave, latexLabel = lm.label_decomposition(chemLabel, scalar_output=True)
#     dict_label = f'{chemLabel}/H1_4861A'
#     theoRatio = H1.getEmissivity(temp, den, wave=wave) / H1.getEmissivity(temp, den, wave=4861)
#     theoEmis_dict[dict_label] = theoRatio
#
# for i, obj in enumerate(objList):
#
#     # Data location
#     cube_address = fitsFolder/fileList[i]
#     objFolder = resultsFolder/obj
#     db_address = objFolder / f'{obj}_database.fits'
#     maskFits_address = objFolder/f'{obj}_masks.fits'
#
#     parameter_fits = objFolder/'gauss_flux.fits'
#     Hbeta_flux = fits.getdata(parameter_fits, 'H1_4861A')
#     hdr_plot = fits.getheader(parameter_fits, 'H1_4861A')
#
#     max_values = [5.0, 5.0, 0.05, 0.05, 0.06, 0.06]
#
#     # Loop through the HI lines
#     HI_lines = list(label_Conver.keys())
#     for j, chemLabel in enumerate(HI_lines):
#
#         line_flux = fits.getdata(parameter_fits, chemLabel)
#
#         ratio_key = f'{chemLabel}/H1_4861A'
#         coeff_im = line_flux/Hbeta_flux
#
#
#         divnorm = colors.TwoSlopeNorm(vmin=0.0,
#                                       vcenter=theoEmis_dict[ratio_key],
#                                       vmax=max_values[j])
#
#         fig = plt.figure(figsize=(10, 10))
#         ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y'))
#         im = ax.imshow(coeff_im, cmap='RdBu', norm=divnorm)
#         cbar = fig.colorbar(im, ax=ax)
#         cbar.set_label(f'Line ratio, white theoretical value ({theoEmis_dict[ratio_key]:.3f})', rotation=270, labelpad=50, fontsize=15)
#         ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_Conver[chemLabel], latex_Conver['H1_4861A'])
#         ax.update({'title': r'Galaxy {}: {}'.format(obj, ratio_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
#         ax.set_xlim(95, 205)
#         ax.set_ylim(75, 225)
#         # plt.show()
#         plt.savefig(objFolder/'extinction'/f'map_{obj}_{chemLabel}_Vabs.png')

# ------ [OIII] line ratios

for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    db_address = objFolder / f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    parameter_fits = objFolder/'intg_flux.fits'
    O3_4959A = fits.getdata(parameter_fits, 'O3_4959A')
    O3_5007A = fits.getdata(parameter_fits, 'O3_5007A')
    hdr_plot = fits.getheader(parameter_fits, 'O3_5007A')

    ion_array, wave_array, latex_array = lime.label_decomposition(['O3_4959A', 'O3_5007A'])

    # coeff_im = O3_5007A/O3_4959A
    #
    # divnorm = colors.TwoSlopeNorm(vmin=2.0,
    #                               vcenter=2.984,
    #                               vmax=4.0)
    # cbar_label = f'Line ratio, theoretical value ({2.984}) white'

    coeff_im = (O3_5007A/O3_4959A - 2.984) * 100

    divnorm = colors.TwoSlopeNorm(vmin=-75.0,
                                  vcenter=0,
                                  vmax=75.0)

    cbar_label = f'Line ratio discrepancy %'

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y'))
    im = ax.imshow(coeff_im, cmap='RdBu', norm=divnorm)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label, rotation=270, labelpad=50, fontsize=15)
    ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_array[1].replace('$', ''), latex_array[0].replace('$', ''))
    ax.update({'title': r'Galaxy {}: {}'.format(obj, ratio_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    ax.set_xlim(95, 205)
    ax.set_ylim(75, 225)
    # plt.show()
    plt.savefig(objFolder/'line_ratios'/f'map_{obj}_OIII_ratio.png')





# # ------ Radial velocities coefficients
# param = 'v_r'
# for i, obj in enumerate(objList):
#
#     # Data location
#     cube_address = fitsFolder/fileList[i]
#     objFolder = resultsFolder/obj
#     db_address = objFolder / f'{obj}_database.fits'
#     maskFits_address = objFolder/f'{obj}_masks.fits'
#
#     parameter_fits = objFolder/f'{param}.fits'
#     hdr_plot = fits.getheader(parameter_fits, 'H1_6563A')
#     vel_Halpha = fits.getdata(parameter_fits, 'H1_6563A')
#     vel_mean = np.nanmean(vel_Halpha)
#     vel_std = np.nanstd(vel_Halpha)
#
#     # Loop through the HI lines
#     HI_lines = list(dinamicLines.keys())
#     for j, chemLabel in enumerate(HI_lines):
#
#         vel_line = fits.getdata(parameter_fits, chemLabel)
#         v_rel = vel_line
#
#         param_min, param_max = np.nanmin(v_rel), np.nanmax(v_rel)
#
#         divnorm = colors.TwoSlopeNorm(vcenter=0.0, vmin=param_min, vmax=param_max)
#
#         fig = plt.figure(figsize=(10, 10))
#         ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y'))
#         im = ax.imshow(v_rel, cmap='RdBu_r', norm=divnorm)
#
#         cbar = fig.colorbar(im)
#         label_bar = latex_Conver[param]
#         cbar.set_label(label_bar, rotation=270, labelpad=50, fontsize=20)
#
#         latex_label = dinamicLines[chemLabel]
#         title_label = f'$v_{{r}}$ {latex_label} -' + r'$\overline{v_{H\alpha}}$' + r' $({:.0f} \pm {:.0f}\,km/s)$'.format(vel_mean, vel_std)
#
#         ax.update({'title': title_label, 'xlabel': r'RA', 'ylabel': r'DEC'})
#         ax.set_xlim(95, 205)
#         ax.set_ylim(75, 225)
#         plt.show()

# # ------ Absolute velocities
# param = 'center'
# c_KMpS = 299792.458
# for i, obj in enumerate(objList):
#
#     # Data location
#     cube_address = fitsFolder / fileList[i]
#     objFolder = resultsFolder / obj
#     db_address = objFolder / f'{obj}_database.fits'
#     maskFits_address = objFolder / f'{obj}_masks.fits'
#
#     parameter_fits = objFolder / f'{param}.fits'
#     hdr_plot = fits.getheader(parameter_fits, 'H1_6563A')
#
#     # Loop through the HI lines
#     HI_lines = list(dinamicLines.keys())
#     for j, chemLabel in enumerate(HI_lines):
#         ion, wave_line, latex_line = lime.label_decomposition(chemLabel, scalar_output=True)
#         center_line_im = fits.getdata(parameter_fits, chemLabel)
#         v_abs = c_KMpS * (1 - wave_line/center_line_im)
#
#         param_mean = np.nanmean(v_abs)
#
#         divnorm = colors.TwoSlopeNorm(vcenter=param_mean, vmin=1375, vmax=1475)
#
#         fig = plt.figure(figsize=(10, 10))
#         ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y'))
#         im = ax.imshow(v_abs, cmap='RdBu_r', norm=divnorm)
#
#         cbar = fig.colorbar(im)
#         label_bar = f'{latex_Conver["v_r"]} (absolute)'
#         cbar.set_label(label_bar, rotation=270, labelpad=50, fontsize=20)
#
#         latex_label = dinamicLines[chemLabel]
#         title_label = f'$v_{{r}}$ {latex_label}' + r' $({:.0f} \pm {:.0f}\,km/s)$'.format(np.nanmean(v_abs), np.nanstd(v_abs))
#
#         ax.update({'title': title_label, 'xlabel': r'RA', 'ylabel': r'DEC'})
#         ax.set_xlim(95, 205)
#         ax.set_ylim(75, 225)
#         # plt.show()
#         plt.savefig(objFolder/'kinematics'/f'map_{obj}_{chemLabel}_Vabs.png')
#
# # ------ Relative velocities wrt Halpha
# param = 'center'
# c_KMpS = 299792.458
# min_bar = [-40, -40, -40, -40, -40, -40, -40, -40]
# max_bar = [50, 50, 50, 50, 50, 50, 50, 50]
# for i, obj in enumerate(objList):
#
#     # Data location
#     cube_address = fitsFolder / fileList[i]
#     objFolder = resultsFolder / obj
#     db_address = objFolder / f'{obj}_database.fits'
#     maskFits_address = objFolder / f'{obj}_masks.fits'
#
#     parameter_fits = objFolder / f'{param}.fits'
#     hdr_plot = fits.getheader(parameter_fits, 'H1_6563A')
#
#     center_Halpha = fits.getdata(parameter_fits, 'H1_6563A')
#     v_abs_Halpha = c_KMpS * (1 - 6563/center_Halpha)
#     v_halpha_mean, v_Halpha_std = np.nanmean(v_abs_Halpha), np.nanstd(v_abs_Halpha)
#
#     # Loop through the HI lines
#     HI_lines = list(dinamicLines.keys())
#     for j, chemLabel in enumerate(HI_lines):
#         ion, wave_line, latex_line = lime.label_decomposition(chemLabel, scalar_output=True)
#         center_line_im = fits.getdata(parameter_fits, chemLabel)
#         v_abs_line = c_KMpS * (1 - wave_line / center_line_im)
#         v_abs_line_std = np.nanstd(v_abs_line)
#
#         v_rel = v_abs_line - v_halpha_mean
#
#         v_mean = np.nanmean(v_abs_line)
#         v_rel_mean = v_mean - v_halpha_mean
#         v_err = np.sqrt(np.power(v_Halpha_std, 2) + np.power(v_abs_line_std,2))
#
#         print(f'{chemLabel} vel_mean = {v_mean:.1f}; Halpha_mean {v_halpha_mean:.1f}; dif {v_mean-v_halpha_mean:.1f}')
#         # print(v_rel-v_halpha_mean)
#
#         print(v_mean-v_halpha_mean, min_bar[j], max_bar[j])
#         divnorm = colors.TwoSlopeNorm(vcenter=v_mean-v_halpha_mean, vmin=min_bar[j], vmax=max_bar[j])
#
#         fig = plt.figure(figsize=(10, 10))
#         ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y'))
#         im = ax.imshow(v_rel, cmap='RdBu_r', norm=divnorm)
#
#         cbar = fig.colorbar(im)
#         label_bar = f'{latex_Conver["v_r"]}' r'(relative to $H\alpha$)'
#         cbar.set_label(label_bar, rotation=270, labelpad=50, fontsize=20)
#         latex_label = dinamicLines[chemLabel]
#         title_label = f'$v_{{r}}$ {latex_label} -' + r'$\overline{v_{H\alpha}}$' + r'= ${:.0f} \pm {:.0f}\,km/s$'.format(v_rel_mean, v_err)
#         ax.update({'title': title_label, 'xlabel': r'RA', 'ylabel': r'DEC'})
#         ax.set_xlim(95, 205)
#         ax.set_ylim(75, 225)
#         # plt.show()
#         plt.savefig(objFolder/'kinematics'/f'map_{obj}_{chemLabel}_Vrel.png')