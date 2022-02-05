import numpy as np
import pandas as pd
import time
import lime
from scipy.interpolate import interp1d
from pathlib import Path
from astropy.io import fits

from src.specsiser.physical_model.extinction_model import ExtinctionModel
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_fado_cube
from progressbar import progressbar

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_abs_list']
fitsFolder_abs = Path(obsData['data_location']['fits_abs_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
norm_flux = obsData['sample_data']['norm_flux_abs']
noise_region = obsData['sample_data']['noiseRegion_abs_array']

dict_errs = {}
dict_nan_values = {}

verbose = False

for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder_abs/fileList[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_addresss = objFolder/f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Output data:
    fitsLog_addresss = objFolder/f'{obj}_linesLog_abs.fits'
    hdul_lineslog = fits.HDUList()

    # Load data
    wave, cube, header = import_fado_cube(cube_address)

    # for idx_region in [0, 1, 2, 3, 4, 5]:
    for idx_region in [0, 1, 2, 3, 4, 5]:

        # Voxel mask
        region_label = f'mask_{idx_region}'
        region_mask = fits.getdata(maskFits_address, region_label, ver=1)
        region_mask = region_mask.astype(bool)
        idcs_voxels = np.argwhere(region_mask)

        # Lines mask
        mask_address = dataFolder/f'{obj}_region_mask_abs.txt'
        mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)

        print(f'\n- Treating {region_label} with {idcs_voxels.shape[0]} pixels')

        # Loop through voxels
        n_lines = 0
        n_voxels = idcs_voxels.shape[0]

        for idx_voxel in np.arange(n_voxels):

            idx_j, idx_i = idcs_voxels[idx_voxel]
            voxel_dict = {}

            local_mask = voxelFolder/f'{idx_j}-{idx_i}_mask_{obj}.txt'
            local_lineslog = voxelFolder/f'{idx_j}-{idx_i}_lineslog_{obj}.txt'
            grid_address_i = voxelFolder/f'{idx_j}-{idx_i}_LineGrid_{obj}.png'

            # Voxel data
            flux_voxel = cube[:, idx_j, idx_i]
            voxel = lime.Spectrum(wave, flux_voxel, norm_flux=norm_flux, crop_waves=(4515, 9500))

            # lime.MaskInspector(mask_address, mask_df, wave, flux_voxel, norm_flux=norm_flux)

            # Locate the line fluxes
            peaks_table, matched_DF = voxel.match_line_mask(mask_df, noise_region, line_type='absorption',
                                                            width_tol=6,
                                                            width_mode='fix')

            if verbose:
                voxel.plot_spectrum(spec_label=f'{obj} voxel {idx_j}-{idx_i}',
                                    peaks_table=peaks_table, matched_DF=matched_DF)

            # Fit and check the regions
            idcsObsLines = (matched_DF.observation == 'detected')
            obsLines = matched_DF.loc[idcsObsLines].index.values
            for j, lineLabel in enumerate(obsLines):
                wave_regions = matched_DF.loc[lineLabel, 'w1':'w6'].values
                try:
                    voxel.fit_from_wavelengths(lineLabel, wave_regions, emission=False, adjacent_cont=False)
                except ValueError as e:
                    err_value = 'NAN values' if 'NaN' in str(e) else 'valueError'
                    err_label = f'ER_{lineLabel[lineLabel.find("_")+1:]}'
                    voxel_dict[err_label] = err_value
                    dict_errs[f'{idx_j}-{idx_i}_{lineLabel}'] = e
                    print(f'--- Line measuring failure at {lineLabel} ({err_value}), {idx_j}-{idx_i}')

            # Spectrum data
            n_lines += len(voxel.log.index)
            voxel_dict['N_Lines'] = len(voxel.log.index)

            # Converting linesLog to fits
            linesHDU = lime.io.lineslog_to_HDU(voxel.log, ext_name=f'{idx_j}-{idx_i}_linelog', header_dict=voxel_dict)

            # Save spectrum data:
            if linesHDU is not None:
                hdul_lineslog.append(linesHDU)

        # Store the drive
        hdul_lineslog.writeto(fitsLog_addresss, overwrite=True, output_verify='fix')
        end = time.time()
