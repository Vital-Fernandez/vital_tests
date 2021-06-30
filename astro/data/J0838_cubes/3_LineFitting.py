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
from astro.data.J0838_cubes.common_methods import import_fits_data

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

verbose = False

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        objFolder = resultsFolder
        cube_address_i = fitsFolder/fileList[i]
        mask_address = dataFolder/obsConf['data_location']['mask_global']
        db_address = objFolder / f'{obj}_database.fits'

        # Output data
        voxelFolder = resultsFolder/obj
        color = 'blue' if 'blue' in obj else 'red'

        # Load the data
        wave, data, header = import_fits_data(cube_address_i, crval3=obsConf[obj]['CRVAL3'], frame_idx=0)
        mask_global_DF = sr.lineslogFile_to_DF(mask_address)

        # Loop throught the line regions
        for idx_region in [0, 1, 2]:

            # Voxel mask
            region_label = f'region_{idx_region}'
            region_mask = fits.getdata(db_address, region_label, ver=1)
            region_mask = region_mask.astype(bool)
            idcs_voxels = np.argwhere(region_mask)
            n_voxels = idcs_voxels.shape[0]

            # Lines mask
            user_conf = obsConf[f'region{idx_region}_line_fitting']

            print(f'\n - Treating {region_label} with {n_voxels} pixels')
            for idx_voxel, idx_pair in enumerate(idcs_voxels):

                idx_j, idx_i = idx_pair

                local_mask = voxelFolder/f'{idx_j}-{idx_i}_mask_{color}.txt'
                local_lineslog = voxelFolder/f'{idx_j}-{idx_i}_lineslog_{color}.txt'
                grid_address_i = voxelFolder/f'{idx_j}-{idx_i}_LineGrid_{color}.png'
                pdfTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable_{color}'
                txtTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable_{color}.txt'

                flux_voxel = data[:, idx_j, idx_i]

                lm = sr.LineMesurer(wave, flux_voxel, redshift=z_list[i], normFlux=norm_flux)
                if verbose:
                    lm.plot_spectrum()

                # Identify the emission lines
                norm_spec = lm.continuum_remover(obsConf[obj]['noiseRegion_array'])
                obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=obsConf[obj]['noiseRegion_array'], intLineThreshold=2.5)
                maskLinesDF = lm.match_lines(obsLinesTable, mask_global_DF, tol=10, find_line_borders=False)

                if verbose:
                    lm.plot_spectrum(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF, specLabel=f'{obj} voxel {idx_j}-{idx_i}')
                    lm.plot_line_mask_selection(maskLinesDF, local_mask, logscale=False)

                # Reset and measure the lines
                lm = sr.LineMesurer(wave, flux_voxel, redshift=z_list[i], normFlux=norm_flux)
                obsLines = maskLinesDF.index.values
                for j, lineLabel in enumerate(obsLines):

                    wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
                    lm.fit_from_wavelengths(lineLabel, wave_regions, user_conf=user_conf)
                    lm.print_results(show_plot=True, show_fit_report=True, log_scale=False)

