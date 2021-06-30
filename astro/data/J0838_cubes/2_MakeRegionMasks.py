import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, red_corr_HalphaHbeta_ratio, store_frame_to_fits
from src.specsiser.print.plot import STANDARD_PLOT
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
import time

# # Declare data and files location
# obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
# objList = obsData['data_location']['object_list']
# fileList = obsData['data_location']['file_list']
# fitsFolder = Path(obsData['data_location']['fits_folder'])
# dataFolder = Path(obsData['data_location']['data_folder'])
# resultsFolder = Path(obsData['data_location']['results_folder'])
#
# z_objs = obsData['sample_data']['z_array']
# pertil_array = obsData['sample_data']['percentil_array']
# noise_region = obsData['sample_data']['noiseRegion_array']
# norm_flux = obsData['sample_data']['norm_flux']
#
# dict_errs = {}
# dict_nan_values = {}
#
# ref_flux_line = 'S3_6312A'

# Declare data and files location
obsConf = sr.loadConfData('J0838_cubes.ini')
fitsFolder = Path(obsConf['data_location']['fits_folder'])
dataFolder = Path(obsConf['data_location']['data_folder'])
resultsFolder = Path(obsConf['data_location']['results_folder'])

fileList = obsConf['data_location']['file_list']
objList = obsConf['data_location']['object_list']
z_list = obsConf['sample_data']['z_array']
norm_flux = obsConf['sample_data']['norm_flux']
percentil_array = obsConf['sample_data']['percentil_array']

ref_flux_line = 'O3_5007A'

# Plot set up
labelsDict = {'xlabel': r'RA',
              'ylabel': r'DEC'}
defaultConf = STANDARD_PLOT.copy()
defaultConf.update(labelsDict)
rcParams.update({})

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        objFolder = resultsFolder
        cube_address_i = fitsFolder/fileList[i]
        mask_address = dataFolder/obsConf['data_location']['mask_global']
        db_address = objFolder / f'{obj}_database.fits'

        # Load the data
        wave, data, header = sr.import_fits_data(cube_address_i, instrument='fits-cube', frame_idx=0)
        mask_global_DF = sr.lineslogFile_to_DF(mask_address)

        # Declare voxels to analyse
        flux5007_image = fits.getdata(db_address, f'{ref_flux_line}_flux', ver=1)
        flux5007_levels = np.nanpercentile(flux5007_image, percentil_array)
        ion, wavelength, latexLabel = sr.label_decomposition(ref_flux_line, scalar_output=True)

        for idx_level, flux_level in enumerate(flux5007_levels):

            if idx_level + 1 < len(flux5007_levels):

                if idx_level == 0:
                    maFlux_image = np.ma.masked_where((flux5007_image >= flux5007_levels[idx_level+1]),
                                                      flux5007_image)
                else:
                    maFlux_image = np.ma.masked_where((flux5007_image <= flux5007_levels[idx_level]) &
                                                      (flux5007_image >= flux5007_levels[idx_level+1]),
                                                      flux5007_image)

                idcs_voxels = np.argwhere(maFlux_image.mask)
                print(f'region {idx_level} ({idcs_voxels.shape[0]} pixels)')

                # if verbose:
                fig = plt.figure(figsize=(12, 8))
                ax = fig.add_subplot()
                ax.update({'title': r'{} galaxy, {} flux'.format(obj, latexLabel), 'xlabel': r'RA', 'ylabel': r'DEC'})
                im = ax.imshow(maFlux_image, cmap=cm.gray)
                plt.show()

                mask_name = f'region_{idx_level}'
                mask_hdu = fits.ImageHDU(name=mask_name, data=maFlux_image.mask.astype(int), ver=1)
                store_frame_to_fits(db_address, fits_hdu=mask_hdu, ext_name=mask_name)