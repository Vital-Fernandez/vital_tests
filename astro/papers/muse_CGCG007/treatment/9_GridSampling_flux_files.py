import numpy as np
import pandas as pd
import lime
from src.specsiser import flambda_calc

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import chemical_lines_indexing
from astro.papers.muse_CGCG007.muse_CGCG007_methods import deredd_fluxes
from astropy.io import fits

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

R_v = obsData['Extinction']['R_v']
red_law = obsData['Extinction']['red_law']

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = f'{objFolder}/chemistry'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    obsLog_addresss = objFolder / f'{obj}_linesLog.fits'
    absLog_addresss = objFolder / f'{obj}_linesLog_abs.fits'
    extinction_map = objFolder/f'{obj}_HI_extinction.fits'

    # Extinction data
    cHbeta_map = fits.getdata(extinction_map, extname='cHbeta')
    errcHbeta_map = fits.getdata(extinction_map, extname='cHbeta_err')

    # Loop throught the line regions
    for idx_region in [0, 1, 2, 3]:
    # for idx_region in [0, 1, 2, 3, 4, 5]:

        # Voxel mask
        region_label = f'mask_{idx_region}'
        region_mask = fits.getdata(maskFits_address, region_label, ver=1)
        region_mask = region_mask.astype(bool)
        idcs_voxels = np.argwhere(region_mask)

        # Simulation time estimation
        n_voxels = idcs_voxels.shape[0]
        print(f'\nTreating {region_label} consisting of {n_voxels}')

        # Output data
        output_inten_files = objFolder / f'region_{idx_region}_gridSampling_intensities.txt'
        output_err_files = objFolder / f'region_{idx_region}_gridSampling_errors.txt'

        # Container for the line intensities
        region_lines = obsData['grid_sampling'][f'region_{idx_region}_line_list']
        intenDF = pd.DataFrame(columns=['mask'] + region_lines)
        errDF = pd.DataFrame(columns=['mask'] + region_lines)

        # Loop through the region voxels
        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            idx_j, idx_i = idx_pair
            ext_ref = f'{idx_j}-{idx_i}_linelog'

            # Load voxel data:
            ext_ref = f'{idx_j}-{idx_i}_linelog'
            obs_log = lime.load_lines_log(obsLog_addresss, ext_ref)
            abs_log = lime.load_lines_log(absLog_addresss, ext_ref)
            print(f'\nTreating voxel {idx_j}-{idx_i}: ({idx_voxel}/{n_voxels})')

            # Correct for extinction
            if 'H1_4861A' in obs_log.index:

                # Get line fluxes:
                linesDF = chemical_lines_indexing(region_lines, obs_log, abs_log, obsData)
                linesDF['obsInt'], linesDF['obsIntErr'] = np.nan, np.nan

                # Input lines
                ion_array, wave_array, latex_array = lime.label_decomposition(linesDF.index.values)

                # Extinction curve
                cHbeta, cHbeta_err = cHbeta_map[idx_j, idx_i], errcHbeta_map[idx_j, idx_i]
                cHbeta = cHbeta if cHbeta > 0 else 0.0
                f_lambda = flambda_calc(wave_array, R_v, red_law)

                # Add values to array
                lineFluxes, lineFluxErr = linesDF['obsFlux'].values, linesDF['obsFluxErr'].values
                obsInt, obsIntErr = deredd_fluxes(lineFluxes, lineFluxErr, cHbeta, cHbeta_err, f_lambda)
                linesDF.loc[:, 'obsInt':'obsIntErr'] = np.array([obsInt, obsIntErr]).T

                # Move to the dataframe
                for line_column in intenDF.columns:
                    if line_column in linesDF.index:
                        intenDF.loc[ext_ref, line_column] = linesDF.loc[line_column, 'obsInt']
                        errDF.loc[ext_ref, line_column] = linesDF.loc[line_column, 'obsIntErr']
                intenDF.loc[ext_ref, 'mask'] = idx_region
                errDF.loc[ext_ref, 'mask'] = idx_region

                # Special voxels treatment
                if idx_region in [0, 1]:
                    intenDF.loc[ext_ref, 'O3_5007A'] = 3 * intenDF.loc[ext_ref, 'O3_4959A']
                    errDF.loc[ext_ref, 'O3_5007A'] = 3 * errDF.loc[ext_ref, 'O3_4959A']

        # Store the drive
        lime.save_line_log(intenDF, output_inten_files)
        lime.save_line_log(errDF, output_err_files)


