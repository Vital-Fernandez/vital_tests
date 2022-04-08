import numpy as np
import pandas as pd
import time
import lime
from scipy.interpolate import interp1d
from pathlib import Path
from astropy.io import fits

from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits
from progressbar import progressbar

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
thres_array = obsData['sample_data']['detect_lim_array']

dict_errs = {}
dict_nan_values = {}

verbose = False

for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_addresss = objFolder/f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Output data:
    fitsLog_addresss = objFolder/f'{obj}_linesLog.fits'
    hdul_lineslog = fits.HDUList()

    # Load data
    wave, cube, header = import_muse_fits(cube_address)

    # Loop throught the line regions

    for idx_region in [0]:

        # Voxel mask
        start = time.time()
        n_lines = 0
        region_label = f'mask_{idx_region}'
        region_mask = fits.getdata(maskFits_address, region_label, ver=1)
        region_mask = region_mask.astype(bool)
        idcs_voxels = np.argwhere(region_mask)

        # Lines mask
        mask_address = dataFolder/f'{obj}_region{idx_region}_mask.txt'
        mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)
        user_conf = obsData[f'region{idx_region}_line_fitting']

        print(f'\n- Treating {region_label} with {idcs_voxels.shape[0]} pixels')

        # Loop through voxels
        n_voxels = idcs_voxels.shape[0]

        # for idx_voxel in progressbar(np.arange(n_voxels), redirect_stdout=True):
        for idx_voxel in np.arange(n_voxels):

            idx_j, idx_i = idcs_voxels[idx_voxel]
            voxel_dict = {}

            grid_address_i = voxelFolder/f'{idx_j}-{idx_i}_LineGrid_{obj}.png'
            pdfTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable'
            txtTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable.txt'

            flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
            flux_err = np.sqrt(cube[:, idx_j, idx_i].var.data) * norm_flux

            voxel = lime.Spectrum(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], norm_flux=norm_flux)

            if verbose:
                voxel.plot_spectrum(spec_label=f'{obj} voxel {idx_j}-{idx_i}', log_scale=True)

            # Security check for pixels with nan values:
            idcs_nan = np.isnan(voxel.flux)
            flux_interpolated = None

            # Nans in flux
            idcs_nan = np.isnan(voxel.flux)
            if idcs_nan.any():
                Interpolation = interp1d(voxel.wave[~idcs_nan], voxel.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
                voxel.flux = Interpolation(voxel.wave)
                if verbose:
                    voxel.plot_spectrum(comp_array=Interpolation(voxel.wave))

            # Nans in errflux
            idcs_nan = np.isnan(voxel.err_flux)
            if idcs_nan.any():
                Interpolation = interp1d(voxel.wave[~idcs_nan], voxel.err_flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
                voxel.err_flux = Interpolation(voxel.wave)

            peaks_table, matched_masks_DF = voxel.match_line_mask(mask_df, noise_region, detect_threshold=thres_array[idx_region],
                                                                  width_mode='fixed')

            # Index of detected lines
            idcsObsLines = matched_masks_DF.index

            if verbose:
                voxel.plot_spectrum(peaks_table=peaks_table, match_log=matched_masks_DF,
                                    spec_label=f'{obj} voxel {idx_j}-{idx_i}', log_scale=True)

            # Fit and check the regions
            obsLines = matched_masks_DF.loc[idcsObsLines].index.values
            for j, lineLabel in enumerate(obsLines):
                # if lineLabel == 'He1_5876A':
                    wave_regions = matched_masks_DF.loc[lineLabel, 'w1':'w6'].values
                    voxel.fit_from_wavelengths(lineLabel, wave_regions, user_cfg=user_conf)
                    voxel.display_results()
                # try:
                #     voxel.fit_from_wavelengths(lineLabel, wave_regions, user_cfg=user_conf)
                #
                # except ValueError as e:
                #     err_value = 'NAN values' if 'NaN' in str(e) else 'valueError'
                #     err_label = f'ER_{lineLabel[lineLabel.find("_")+1:]}'
                #     voxel_dict[err_label] = err_value
                #     dict_errs[f'{idx_j}-{idx_i}_{lineLabel}'] = e
                #     print(f'--- Line measuring failure at {lineLabel} ({err_value}), {idx_j}-{idx_i}')

            # Spectrum data
            n_lines += len(voxel.log.index)
            voxel_dict['N_Lines'] = len(voxel.log.index)
            voxel_dict['N_nan'] = idcs_nan.sum()

            # Converting linesLog to fits
            linesHDU = lime.io.log_to_HDU(voxel.log, ext_name=f'{idx_j}-{idx_i}_linelog', header_dict=voxel_dict)

            # Save spectrum data:
            if linesHDU is not None:
                hdul_lineslog.append(linesHDU)

        # Store the drive
        # hdul_lineslog.writeto(fitsLog_addresss, overwrite=True, output_verify='fix')
        end = time.time()
        print(f'- Execution time {(end - start):.2f} seconds, {(end - start) / 60:.1f} min, for {n_lines} lines, voxels {n_voxels}')
        # print(
        #     f'- Execution time {(end - start):.2f} seconds, {(end - start) / 60:.1f} min, for {n_lines} lines, voxels {n_voxels},'
        #     f' errors {len(dict_errs.keys())}')

# Show summary
print('-- Error summary')
for voxel_fail, error in dict_errs.items():
    print(voxel_fail)
