import numpy as np
import pandas as pd
import time
import lime
from src.specsiser import flambda_calc

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import HII_chemistry_label_conversion as HII_conv, epm_HII_CHI_mistry_orig
from astro.papers.muse_CGCG007.muse_CGCG007_methods import deredd_fluxes
from astropy.io import fits

pd.set_option('display.max_columns', None)

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

list_errors = []

HII_chim_folder = Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v5.22/')

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = f'{objFolder}/chemistry'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    obsLog_addresss = objFolder / f'{obj}_linesLog.fits'
    extinction_map = objFolder/f'{obj}_HI_extinction.fits'

    # Extinction data
    cHbeta_map = fits.getdata(extinction_map, extname='cHbeta')
    errcHbeta_map = fits.getdata(extinction_map, extname='cHbeta_err')

    # Output data
    HIIchimistry_fits = objFolder / f'{obj}_HIIchimistry.fits'
    hdul_lineslog = fits.HDUList()

    # Loop throught the line regions
    start = time.time()
    for idx_region in [0, 1]:

        # Voxel mask
        region_label = f'mask_{idx_region}'
        region_mask = fits.getdata(maskFits_address, region_label, ver=1)
        region_mask = region_mask.astype(bool)
        idcs_voxels = np.argwhere(region_mask)

        # Simulation time estimation
        n_voxels = idcs_voxels.shape[0]
        print(f'\nTreating {region_label} consisting of {n_voxels}')

        # Store in HII-chims-try format
        HII_CHI_mistry_DF = pd.DataFrame()

        # Loop through the region voxels
        n_voxels = idcs_voxels.shape[0]
        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            idx_j, idx_i = idx_pair
            print(f'\nTreating voxel {idx_j}-{idx_i}: ({idx_voxel}/{n_voxels})')

            # Load voxel data:
            ext_ref = f'{idx_j}-{idx_i}_linelog'
            obs_log = lime.load_lines_log(obsLog_addresss, ext_ref)
            cHbeta, cHbeta_err = cHbeta_map[idx_j, idx_i], errcHbeta_map[idx_j, idx_i]

            # Prepare data for HII-CHI-mistry
            idcs_inputLines = obs_log.index.isin(HII_conv.keys())
            input_lines = obs_log.loc[idcs_inputLines].index.values
            ion, wave, latex = lime.label_decomposition(input_lines)

            # Correct for extinction
            if 'H1_4861A' in obs_log.index:

                flux_Hbeta = obs_log.loc['H1_4861A', 'gauss_flux']
                gauss_fluxes = obs_log.loc[idcs_inputLines].gauss_flux.values/flux_Hbeta
                gauss_errs = obs_log.loc[idcs_inputLines].gauss_err.values/flux_Hbeta
                f_lambda = flambda_calc(wave, R_v, red_law)
                obsInt, obsIntErr = deredd_fluxes(gauss_fluxes, gauss_errs, cHbeta, cHbeta_err, f_lambda)
                obs_log['obsInt'], obs_log['obsIntErr'] = np.nan, np.nan
                obs_log.loc[idcs_inputLines, 'obsInt':'obsIntErr'] = np.array([obsInt, obsIntErr]).T

                print(ext_ref)
                print('Input fluxes', gauss_fluxes)
                print('f_lambda', f_lambda)
                print('cHbeta', cHbeta)
                print('cHbeta_err', cHbeta_err)



                # Convert the line names
                for lineLabel in input_lines:
                    HII_CHI_mistry_DF.loc[ext_ref, HII_conv[lineLabel]] = obs_log.loc[lineLabel, 'obsInt']
                    HII_CHI_mistry_DF.loc[ext_ref, f'e{HII_conv[lineLabel]}'] = obs_log.loc[lineLabel, 'obsIntErr']

                # Region selectiontreatment:
                if idx_region in [0, 1]:
                    HII_CHI_mistry_DF.loc[ext_ref, 'OIII_5007'] = 3 * HII_CHI_mistry_DF.loc[ext_ref, 'OIII_4959']
                    HII_CHI_mistry_DF.loc[ext_ref, 'eOIII_5007'] = 3 * HII_CHI_mistry_DF.loc[ext_ref, 'eOIII_4959']

                # # Run HII-CHI-mistry
                # try:
                #     lineSA = HII_CHI_mistry_DF.to_records(index=False)
                #     outputSA = epm_HII_CHI_mistry(lineSA, output_file='None', n=200, sed=1, inter=1)
                #     linesCol = fits.ColDefs(outputSA)
                #     linesHDU = fits.BinTableHDU.from_columns(linesCol, name=f'{idx_j}-{idx_i}_HIIchimistry')
                #     hdul_lineslog.append(linesHDU)
                # except:
                #     list_errors.append(ext_ref)

        # Store the drive
        HII_CHI_mistry_DF.reset_index(inplace=True)
        HII_CHI_mistry_DF.rename(columns={'index': 'ID'}, inplace=True)
        region_HII_chim_lines_file = HII_chim_folder/f'region_{idx_region}_fluxes.txt'
        with open(region_HII_chim_lines_file, 'wb') as output_file:
            string_DF = HII_CHI_mistry_DF.to_string(index=False)
            output_file.write(string_DF.encode('UTF-8'))

        # Run command example:
        # python HCm_v5.22.py region_0_fluxes.txt 200


    #     hdul_lineslog.writeto(HIIchimistry_fits, overwrite=True, output_verify='fix')
    # end = time.time()
#
# print(f'- Execution time {(end - start)/60:.3f} min')
#
# print(f'These {len(list_errors)} voxels failed:')
# for err_vox in list_errors:
#     print(err_vox)

# 'python HCm v5.2.py input.txt 100'
