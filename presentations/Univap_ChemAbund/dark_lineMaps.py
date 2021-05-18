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

conf_file = Path('/home/vital/PycharmProjects/vital_tests/astro/data/muse/muse_greenpeas.ini')

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

        # Recombination lines
        for chemLabel, plotLabel in label_Conver.items():

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))

            dict_label = f'{plotLabel}/Hbeta'
            flux_image = fits.getdata(linesMaps_fits_address, dict_label, ver=1)

            divnorm = colors.TwoSlopeNorm(vmin=np.nanmin(flux_image),
                                          vcenter=theoEmis_dict[dict_label],
                                          vmax=np.nanmax(flux_image))

            im = ax.imshow(flux6563_image, cmap=halpha_cmap,
                           norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))
            im2 = ax.imshow(flux_image, cmap='RdBu', norm=divnorm)

            cbar = fig.colorbar(im2, ax=ax)
            cbar.set_label('Line ratio (white theoretical value)', rotation=270, labelpad=50, fontsize=15)

            ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_Conver[chemLabel], latex_Conver['H1_4861A'])
            ax.update({'title': r'Galaxy {}: {}'.format(obj, ratio_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
            plt.tight_layout()
            plt.show()

        # Recombination lines
        for dinLabel, latex_label in dinamicLines.items():
            for param in ('v_r', 'sigma_vel'):

                fig = plt.figure(figsize=(10, 10))
                ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))

                dict_label = f'{param}_{dinLabel}'
                param_image = fits.getdata(linesMaps_fits_address, dict_label, ver=1)

                im = ax.imshow(flux6563_image, cmap=halpha_cmap,
                               norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))

                if param == 'v_r':
                    # param_min, param_max = np.nanmin(param_image), np.nanmax(param_image)
                    # divnorm = colors.TwoSlopeNorm(vmin=param_min,
                    #                               vcenter=0.0,
                    #                               vmax=param_max)
                    # im2 = ax.imshow(param_image, cmap='RdBu', norm=divnorm)
                    im2 = ax.imshow(param_image)

                else:
                    im2 = ax.imshow(param_image)

                cbar = fig.colorbar(im2)
                label_bar = latex_Conver[param]
                cbar.set_label(label_bar, rotation=270, labelpad=50, fontsize=20)

                ax.update({'title': r'Galaxy {}: {}'.format(obj, latex_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
                plt.tight_layout()
                # plt.savefig(resultsFolder/obj/f'map_{obj}_{dinLabel}_{param}.png', bbox_inches='tight')
                plt.show()