import numpy as np
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import background_color, DARK_PLOT, label_Conver, latex_Conver, dinamicLines
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astropy.table import Table
import pyneb as pn
from src.specsiser.print.plot import STANDARD_PLOT

sulfur_lines = {'S2_6716A': r'$[SII]6716\AA$', 'S2_6731A': r'$[SII]6716\AA$'}

# Plot set up
defaultConf = STANDARD_PLOT.copy()
defaultConf['axes.titlesize'] = 20
rcParams.update(defaultConf)

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
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

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        objFolder = resultsFolder/obj
        voxelFolder = resultsFolder/obj/'voxel_data'
        db_address = objFolder / f'{obj}_database.fits'
        fitsLog_address = objFolder / f'{obj}_linesLog.fits'

        # Output file
        linesMaps_fits_address = objFolder/f'{obj}_lineParamMaps.fits'

        # Extinction model
        red_model = sr.ExtinctionModel(Rv=obsData['Extinction']['R_v'], red_curve=obsData['Extinction']['red_law'])

        regions_to_treat = [0, 1, 2, 3]

        # # ----------------------------------------- Generate the image data
        #
        # # Empty containers for the images
        # image_dict = {}
        # for chemLabel, plotLabel in label_Conver.items():
        #     dict_label = f'{plotLabel}/Hbeta'
        #     image_dict[dict_label] = np.full(voxel_grid_size.astype(int), np.nan)
        #
        # for dinLabel, plotLabel in dinamicLines.items():
        #     for param in ('v_r', 'sigma_vel'):
        #         image_dict[f'{param}_{dinLabel}'] = np.full(voxel_grid_size.astype(int), np.nan)
        #
        # for dinLabel, plotLabel in sulfur_lines.items():
        #     for param in ('v_r', 'sigma_vel'):
        #         image_dict[f'{param}_{dinLabel}'] = np.full(voxel_grid_size.astype(int), np.nan)
        #
        # # Open the lines log database
        # with fits.open(fitsLog_address) as hdul:
        #
        #     # Loop throught the line regions
        #     for idx_region in regions_to_treat:
        #         region_label = f'region_{idx_region}'
        #         region_mask = fits.getdata(db_address, region_label, ver=1)
        #         region_mask = region_mask.astype(bool)
        #         idcs_voxels = np.argwhere(region_mask)
        #
        #         # Loop through the region voxels
        #         for idx_voxel, idx_pair in enumerate(idcs_voxels):
        #             idx_j, idx_i = idx_pair
        #             logLabel = f'{idx_j}-{idx_i}_linelog'
        #
        #             # Load lines log data and store it as an image
        #             if logLabel in hdul:
        #                 linesDF = Table(hdul[logLabel].data).to_pandas()
        #                 linesDF.set_index('index', inplace=True)
        #
        #                 if 'H1_4861A' in linesDF.index:
        #                     Hbeta_flux = linesDF.loc['H1_4861A', 'gauss_flux']
        #                     for chemLabel, plotLabel in label_Conver.items():
        #
        #                         if chemLabel in linesDF.index:
        #                             dict_label = f'{plotLabel}/Hbeta'
        #                             lineFlux = linesDF.loc[chemLabel, 'gauss_flux']/Hbeta_flux
        #                             image_dict[dict_label][idx_j, idx_i] = lineFlux# - theoEmis_dict[dict_label]
        #
        #                 for dinLabel, plotLabel in dinamicLines.items():
        #                     if dinLabel in linesDF.index:
        #                         for param in ('v_r', 'sigma_vel'):
        #                             image_dict[f'{param}_{dinLabel}'][idx_j, idx_i] = linesDF.loc[dinLabel, param]
        #
        # # Storing the dictionary as a fits image file
        # new_hdul = fits.HDUList()
        # new_hdul.append(fits.PrimaryHDU())
        #
        # # Second page for the fits file plot configuration
        # hdr_plot = fits.getheader(db_address, extname='PlotConf')
        # hdu_table = fits.BinTableHDU.from_columns(columns=[], header=hdr_plot, name='PlotConf')
        # new_hdul.append(hdu_table)
        #
        # for param_line, param_map in image_dict.items():
        #     new_hdul.append(fits.ImageHDU(name=param_line, data=param_map, ver=1))
        # new_hdul.writeto(linesMaps_fits_address, overwrite=True, output_verify='fix')
        #
        # # Open the lines log database
        # fitsSIILog_address = objFolder / f'{obj}_linesLog_SII.fits'
        # with fits.open(fitsSIILog_address) as hdul:
        #
        #     # Loop throught the line regions
        #     for idx_region in regions_to_treat:
        #         region_label = f'region_{idx_region}'
        #         region_mask = fits.getdata(db_address, region_label, ver=1)
        #         region_mask = region_mask.astype(bool)
        #         idcs_voxels = np.argwhere(region_mask)
        #
        #         # Loop through the region voxels
        #         for idx_voxel, idx_pair in enumerate(idcs_voxels):
        #             idx_j, idx_i = idx_pair
        #             logLabel = f'{idx_j}-{idx_i}_linelog'
        #
        #             # Load lines log data and store it as an image
        #             if logLabel in hdul:
        #                 linesDF = Table(hdul[logLabel].data).to_pandas()
        #                 linesDF.set_index('index', inplace=True)
        #
        #                 for dinLabel, plotLabel in sulfur_lines.items():
        #                     if dinLabel in linesDF.index:
        #                         for param in ('v_r', 'sigma_vel'):
        #                             image_dict[f'{param}_{dinLabel}'][idx_j, idx_i] = linesDF.loc[dinLabel, param]
        #
        # # Storing the dictionary as a fits image file
        # new_hdul = fits.HDUList()
        # new_hdul.append(fits.PrimaryHDU())
        #
        # # Second page for the fits file plot configuration
        # hdr_plot = fits.getheader(db_address, extname='PlotConf')
        # hdu_table = fits.BinTableHDU.from_columns(columns=[], header=hdr_plot, name='PlotConf')
        # new_hdul.append(hdu_table)
        #
        # for param_line, param_map in image_dict.items():
        #     new_hdul.append(fits.ImageHDU(name=param_line, data=param_map, ver=1))
        # new_hdul.writeto(linesMaps_fits_address, overwrite=True, output_verify='fix')


        # ----------------------------------------- Generate the image plots
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
        #
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
        #     ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_Conver[chemLabel], latex_Conver['H1_4861A'])
        #     ax.update({'title': r'Galaxy {}: {}'.format(obj, ratio_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
        #     plt.tight_layout()
        #     plt.show()

        # Line kinematics
        dinamicLines.update(sulfur_lines)
        vr_Halpha = fits.getdata(linesMaps_fits_address, 'v_r_H1_6563A', ver=1)
        Halpha_label = dinamicLines['H1_6563A']
        Halpha_mean, Halpha_std = np.nanmean(vr_Halpha), np.nanstd(vr_Halpha)
        print(Halpha_mean, Halpha_std)
        for dinLabel, latex_label in dinamicLines.items():
            for param in ['v_r']:#('v_r', 'sigma_vel'):

                fig = plt.figure(figsize=(10, 10))
                fig.tight_layout()

                ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))

                dict_label = f'{param}_{dinLabel}'
                param_image = fits.getdata(linesMaps_fits_address, dict_label, ver=1)

                im = ax.imshow(flux6563_image, cmap=halpha_cmap,
                               norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))

                if param == 'v_r':
                    param_image = param_image - Halpha_mean
                    param_min, param_max = np.nanmin(param_image), np.nanmax(param_image)
                    divnorm = colors.TwoSlopeNorm(vcenter=0.0, vmin=-30, vmax=30)
                    im2 = ax.imshow(param_image, cmap='RdBu_r', norm=divnorm)
                    # title_label = f'$v_{{r}}$ {latex_label} - {Halpha_label}'
                    title_label = f'$v_{{r}}$ {latex_label} -' + r'$\overline{v_{H\alpha}}$' + r' $({:.0f} \pm {:.0f}\,km/s)$'.format(Halpha_mean, Halpha_std)

                else:
                    im2 = ax.imshow(param_image)
                    title_label = f'$\sigma_{{int}})$ {Halpha_label}'

                cbar = fig.colorbar(im2)
                label_bar = latex_Conver[param]
                cbar.set_label(label_bar, rotation=270, labelpad=50, fontsize=20)

                ax.update({'title': r'Galaxy {}: {}'.format(obj, title_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
                ax.set_xlim(100, 210)
                ax.set_ylim(90, 240)

                # plt.subplots_adjust(top=0.85)
                plt.savefig(resultsFolder/obj/f'map_{obj}_{dinLabel}_{param}.png')
                # plt.show()