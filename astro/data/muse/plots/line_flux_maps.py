import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, red_corr_HalphaHbeta_ratio, default_linelog_types
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
import time
from astropy.table import Table
import pyneb as pn
from src.specsiser.print.plot import STANDARD_PLOT
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

dict_errs = {}
dict_nan_values = {}

ref_flux_line = 'S3_6312A'
sulfur_bdry = int(obsData['Region_masks']['S3_direct_limit'])
hydrogen_bdry = int(obsData['Region_masks']['H1_direct_limit'])
verbose = False

H1 = pn.RecAtom('H', 1)
temp, den = 10000.0, 100.0

label_Conver = {'H1_6563A': 'Halpha',
               'H1_9229A': 'HPas9',
               'H1_9015A': 'HPas10',
               'H1_8863A': 'HPas11',
               'H1_8750A': 'HPas12'}

latex_Conver = {'H1_6563A': r'H\alpha',
                'H1_4861A': r'H\beta',
                'H1_9229A': r'H_{Pas,\,9}',
                'H1_9015A': r'H_{Pas,\,10}',
                'H1_8863A': r'H_{Pas,\,11}',
                'H1_8750A': r'H_{Pas,\,12}',
                'v_r': r'$v_{r}\,(km/s)$',
                'sigma_vel': r'$\sigma_{int}\,(km/s)$'}

dinamicLines = {'H1_6563A': r'$H\alpha_{Narrow}$',
              'H1_6563A_w1': r'$H\alpha_{Broad\,1}$',
              'H1_6563A_w2': r'$H\alpha_{Broad\,2}$',
              'H1_4861A': r'$H\beta_{Narrow}$',
              'H1_4861A_w1': r'$H\beta_{Broad\,1}$',
              'O3_5007A': r'$[OIII]5007\AA_{Narrow}$',
              'O3_5007A_w1': r'$[OIII]5007\AA_{Broad\,1}$',
              'O3_5007A': r'$[OIII]5007\AA_{Narrow}$',
              'O3_5007A_w1': r'$[OIII]5007\AA_{Broad\,1}$'}

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
        db_addresss = objFolder / f'{obj}_database.fits'
        fitsLog_addresss = objFolder / f'{obj}_linesLog.fits'

        # Extinction model
        red_model = sr.ExtinctionModel(Rv=obsData['Extinction']['R_v'], red_curve=obsData['Extinction']['red_law'])

        # Empty containers for the images
        image_dict = {}
        for chemLabel, plotLabel in label_Conver.items():
            dict_label = f'{plotLabel}/Hbeta'
            image_dict[dict_label] = np.full(voxel_grid_size.astype(int), np.nan)

        kinema_dict = {'v_r': {}, 'sigma_vel': {}}
        for dinLabel, plotLabel in dinamicLines.items():
            kinema_dict['v_r'][dinLabel] = np.full(voxel_grid_size.astype(int), np.nan)
            kinema_dict['sigma_vel'][dinLabel] = np.full(voxel_grid_size.astype(int), np.nan)

        # Open the lines log database
        with fits.open(fitsLog_addresss) as hdul:

            # Loop throught the line regions
            for idx_region in [0, 1, 2]:

                # Voxel mask
                region_label = f'region_{idx_region}'
                region_mask = fits.getdata(db_addresss, region_label, ver=1)
                region_mask = region_mask.astype(bool)
                idcs_voxels = np.argwhere(region_mask)

                for idx_voxel, idx_pair in enumerate(idcs_voxels):

                    idx_j, idx_i = idx_pair
                    logLabel = f'{idx_j}-{idx_i}_linelog'

                    if logLabel in hdul:

                        linesDF = Table(hdul[logLabel].data).to_pandas()
                        linesDF.set_index('index', inplace=True)

                        if 'H1_4861A' in linesDF.index:

                            Hbeta_flux = linesDF.loc['H1_4861A', 'gauss_flux']

                            for chemLabel, plotLabel in label_Conver.items():

                                if chemLabel in linesDF.index:
                                    dict_label = f'{plotLabel}/Hbeta'
                                    lineFlux = linesDF.loc[chemLabel, 'gauss_flux']/Hbeta_flux
                                    image_dict[dict_label][idx_j, idx_i] = lineFlux# - theoEmis_dict[dict_label]

                        for dinLabel, plotLabel in dinamicLines.items():

                            if dinLabel in linesDF.index:

                                for param in ('v_r', 'sigma_vel'):
                                    kinema_dict[param][dinLabel][idx_j, idx_i] = linesDF.loc[dinLabel, param]

        # Generate the image maps
        hdr_plot = fits.getheader(db_addresss, extname='PlotConf')
        flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
        min_level = fits.getval(db_addresss, keyword=f'P8000', extname=f'H1_6563A_flux')

        # # ------------------------- Recombination line ratios
        # for chemLabel, plotLabel in label_Conver.items():
        #
        #     fig = plt.figure(figsize=(10, 10))
        #     ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))
        #
        #     dict_label = f'{plotLabel}/Hbeta'
        #     flux_image = image_dict[dict_label]
        #     print(np.nanmin(flux_image), np.nanmax(flux_image))
        #
        #     divnorm = colors.TwoSlopeNorm(vmin=np.nanmin(flux_image), vcenter=theoEmis_dict[dict_label], vmax=np.nanmax(flux_image))
        #     im = ax.imshow(flux6563_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=min_level, vmin=min_level, base=10))
        #     im2 = ax.imshow(flux_image, cmap='RdBu', norm=divnorm)
        #
        #     cbar = fig.colorbar(im2, ax=ax)
        #     cbar.set_label('Line ratio (white theoretical value)', rotation=270, labelpad=50, fontsize=15)
        #
        #     ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_Conver[chemLabel], latex_Conver['H1_4861A'])
        #     ax.update({'title': r'Galaxy {}: {}'.format(obj, ratio_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
        #     plt.tight_layout()
        #     # plt.savefig(resultsFolder/obj/f'map_{obj}_{parameter}.png', )
        #     plt.show()

        # ------------------------- Line kinematics
        for dinLabel, latex_label in dinamicLines.items():
            for param in ('v_r', 'sigma_vel'):

                fig = plt.figure(figsize=(10, 10))
                ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))

                param_image = kinema_dict[param][dinLabel]
                param_min, param_max = np.nanmin(param_image), np.nanmax(param_image)
                # param_cen = 0 if param == 'v_r' else np.nanmean(param_image)
                # print(dinLabel, param)
                # print(param_min, param_cen, param_max)
                # print()
                # divnorm = colors.TwoSlopeNorm(vmin=param_min, vcenter=param_cen, vmax=param_max)

                im = ax.imshow(flux6563_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=min_level, vmin=min_level, base=10))
                im2 = ax.imshow(param_image)
                cbar = fig.colorbar(im2)
                label_bar = latex_Conver[param]

                # label_bar = r'$\frac{H\alpha}{H\beta} - \left(\frac{H\alpha}{H\beta}\right)_{T_{e}=10000\,K,\,n_{e}=100cm^{-3}}$'
                cbar.set_label(label_bar, rotation=270, labelpad=50, fontsize=20)

                ax.update({'title': r'Galaxy {}: {}'.format(obj, latex_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
                plt.tight_layout()
                # plt.savefig(resultsFolder/obj/f'map_{obj}_{dinLabel}_{param}.png', bbox_inches='tight')
                plt.show()