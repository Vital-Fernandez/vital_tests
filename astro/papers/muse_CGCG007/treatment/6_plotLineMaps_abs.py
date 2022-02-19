import numpy as np

import lime
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS

from astro.papers.muse_CGCG007.muse_CGCG007_methods import label_Conver, abs_target_lines, param_images_abs, abs_target_lines

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
voxel_grid_size = obsData['sample_data']['grid_shape_array']
coordinates_keys_list = obsData['data_location']['wcs_key_list']

# # ------ Absorption line ratios
# for i, obj in enumerate(objList):
#
#     # Data location
#     cube_address = fitsFolder/fileList[i]
#     objFolder = resultsFolder/obj
#     db_address = objFolder / f'{obj}_database.fits'
#     maskFits_address = objFolder/f'{obj}_masks.fits'
#
#     parameter_fits = objFolder/'abs_gauss_flux.fits'
#     Hbeta_flux = fits.getdata(parameter_fits, 'H1_4861A')
#     hdr_plot = fits.getheader(parameter_fits, 'H1_4861A')
#
#     ion_Hbeta, wave_Hbeta, latex_Hbeta = lime.label_decomposition('H1_4861A', scalar_output=True)
#     ion_array, wave_array, latex_array = lime.label_decomposition(abs_target_lines)
#
#     for j, line in enumerate(abs_target_lines):
#
#         line_flux = fits.getdata(parameter_fits, line)
#
#         coeff_im = (line_flux/Hbeta_flux)
#
#         divnorm = colors.TwoSlopeNorm(vmin=-75.0,
#                                       vcenter=0,
#                                       vmax=75.0)
#
#         cbar_label = f'Flux ratio %'
#
#         fig = plt.figure(figsize=(10, 10))
#         ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y'))
#         # im = ax.imshow(coeff_im, cmap='RdBu', norm=divnorm)
#         im = ax.imshow(coeff_im)
#         cbar = fig.colorbar(im, ax=ax)
#         cbar.set_label(cbar_label, rotation=270, labelpad=50, fontsize=15)
#         ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_array[j].replace('$', ''), latex_Hbeta.replace('$', ''))
#         ax.update({'title': r'Galaxy {}: {}'.format(obj, ratio_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
#         ax.set_xlim(95, 205)
#         ax.set_ylim(75, 225)
#         plt.show()
#         plt.savefig(objFolder/'line_ratios'/f'map_abs_ratio_{line}-H1_4861A_ratio.png')

# # ------ Emission/absroption
# for i, obj in enumerate(objList):
#
#     # Data location
#     cube_address = fitsFolder/fileList[i]
#     objFolder = resultsFolder/obj
#     db_address = objFolder / f'{obj}_database.fits'
#     maskFits_address = objFolder/f'{obj}_masks.fits'
#
#     emis_fits = objFolder / f'gauss_flux.fits'
#     abs_fits = objFolder/'abs_gauss_flux.fits'
#     Hbeta_emis = fits.getdata(emis_fits, 'H1_4861A')
#     Hbeta_abs = fits.getdata(abs_fits, 'H1_4861A')
#
#     hdr_plot = fits.getheader(abs_fits, 'H1_4861A')
#
#     for j, line in enumerate(abs_target_lines):
#
#         ion_array, wave_array, latex_array = lime.label_decomposition(line, scalar_output=True)
#
#         cont_emis = fits.getdata(objFolder/'cont.fits', line)
#         cont_abs = fits.getdata(objFolder/'abs_cont.fits', line)
#
#         norm = cont_emis/cont_abs
#
#         line_flux_emis = fits.getdata(emis_fits, line)
#         line_flux_abs = fits.getdata(abs_fits, line) * -1
#
#         # idxY, idxX = 167, 167
#         # norm[idxY, idxX]
#         # line_flux_emis[idxY, idxX]
#         # line_flux_abs[idxY, idxX]
#         # coeff_im[idxY, idxX]
#
#         coeff_im = (line_flux_emis/(line_flux_abs * norm))
#
#         print(np.sum(coeff_im<0))
#
#         fig = plt.figure(figsize=(10, 10))
#         ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y'))
#         im = ax.imshow(coeff_im, norm=colors.LogNorm())
#         cbar = fig.colorbar(im, ax=ax)
#         title = f'{latex_array} ' + r'$\frac{emission}{absorption}$'
#         ax.update({'title': r'Galaxy {}: {}'.format(obj, title), 'xlabel': r'RA', 'ylabel': r'DEC'})
#         ax.set_xlim(95, 205)
#         ax.set_ylim(75, 225)
#         plt.show()
#         # plt.savefig(objFolder/'line_ratios'/f'map_{obj}_OIII_ratio.png')


# ------ Emission/absroption distribution
for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    db_address = objFolder / f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    emis_fits = objFolder / f'gauss_flux.fits'
    abs_fits = objFolder/'abs_gauss_flux.fits'
    Hbeta_emis = fits.getdata(emis_fits, 'H1_4861A')
    Hbeta_abs = fits.getdata(abs_fits, 'H1_4861A')

    hdr_plot = fits.getheader(abs_fits, 'H1_4861A')

    colorNorm = colors.Normalize(0, 6)
    cmap = cm.get_cmap(name=None)

    for j, line in enumerate(abs_target_lines):

        ion_array, wave_array, latex_array = lime.label_decomposition(line, scalar_output=True)

        array_container = []
        data_labels = []
        colors = []
        for idx_region in [0, 1, 2, 3, 4, 5]:

            # Voxel mask
            region_label = f'MASK_{idx_region}'
            region_mask = fits.getdata(maskFits_address, region_label, ver=1)
            region_mask = region_mask.astype(bool)
            # idcs_voxels = np.argwhere(region_mask)
            print(region_label, np.sum(region_mask))

            cont_emis = fits.getdata(objFolder/'cont.fits', line)
            cont_abs = fits.getdata(objFolder/'abs_cont.fits', line)
            line_flux_emis = fits.getdata(emis_fits, line)
            line_flux_abs = fits.getdata(abs_fits, line) * -1

            cont_emis = cont_emis[region_mask]
            cont_abs = cont_abs[region_mask]
            line_flux_emis = line_flux_emis[region_mask]
            line_flux_abs = line_flux_abs[region_mask]

            norm = cont_emis/cont_abs

            coeff_im = (line_flux_emis/(line_flux_abs * norm))
            idcs_non_nan = ~np.isnan(coeff_im) & (coeff_im > 0)

            idcs_neg = coeff_im < 0

            if idcs_neg.size > 0:
                data_array = coeff_im[idcs_non_nan]
                print(f'- {region_label} {np.sum(idcs_neg)} {data_array.shape}')
                array_container.append(data_array)
                data_labels.append(f'{region_label} ({len(data_array)})')
                colors.append(cmap(colorNorm(idx_region)))

        # Computing mean absorptions
        data_total = np.concatenate(array_container)
        mean_abs, std_abs = data_total.mean(), data_total.std()

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.hist(array_container, label=data_labels, bins=30, log=True, stacked=True, color=colors)
        title = f'{latex_array} ' + r'$\frac{emission}{absorption}$'
        ax.axvline(x=mean_abs, color='black', linestyle='--')
        ax.axvspan(mean_abs - std_abs, mean_abs + std_abs, alpha=0.25, color='grey')
        ax.update({'title': r'Galaxy {}'.format(obj), 'xlabel': title, 'ylabel': r'Count'})
        ax.legend()
        plt.show()




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