import numpy as np
import pandas as pd
import time
import lime
from scipy.interpolate import interp1d
from pathlib import Path
from astropy.io import fits
from mpdaf.obj import Cube
from astropy.wcs import WCS


def import_muse_fits(file_address):

    cube = Cube(filename=str(file_address))
    header = cube.data_header

    cube.wave.info()
    dw = header['CD3_3']
    w_min = header['CRVAL3']
    nPixels = header['NAXIS3']
    w_max = w_min + dw * nPixels
    wave = np.linspace(w_min, w_max, nPixels, endpoint=False)

    return wave, cube, header


# Load observation data and configuration
obs_config = lime.load_cfg('../muse_CGCG007.ini')

# Declare data and files location
fits_folder = Path(obs_config['data_location']['fits_folder'])
data_folder = Path(obs_config['data_location']['data_folder'])
results_folder = Path(obs_config['data_location']['results_folder'])
file_list = obs_config['data_location']['file_list'][0]

# Detection parameters
pertil_array = obs_config['sample_data']['percentil_array']
noise_region = obs_config['sample_data']['noiseRegion_array']
thres_array = obs_config['sample_data']['detect_lim_array']

# Observation information
obj_name = obs_config['data_location']['object_list'][0]
redshift = obs_config['sample_data']['z_array'][0]
norm_flux = obs_config['sample_data']['norm_flux']

# Input/outputs paths
cube_address = fits_folder/file_list
mask_address = results_folder/obj_name/f'{obj_name}_masks_SN.fits'
bands_file_0 = data_folder/f'{obj_name}_region{0}_mask.txt'
bands_df_0 = lime.load_log(bands_file_0)
logs_address = results_folder/obj_name/f'{obj_name}_log_lines.fits'

# Load spectrum
print(cube_address)
print(mask_address)

wave, cube, header = import_muse_fits(cube_address)
wcs = WCS(header)

flux = cube.data.data * norm_flux
err = np.sqrt(cube.var.data) * norm_flux

# Declare cube data
obj = lime.Cube(wave, flux, err, redshift, norm_flux)

# Compute the spatial mask
percentiles = (95.50, 97.50, 99.50, 99.90, 99.99)
Halpha_band = bands_df_0.loc['H1_6563A_b', 'w1':'w6'].to_numpy(float)
obj.plot.cube('H1_6563A_b', Halpha_band, line_fg='H1_6563A_b', band_fg=Halpha_band, percentiles_fg=percentiles, wcs=wcs)

obj.spatial_masker('H1_6563A_b', band=Halpha_band, param='SN_line', percentiles=percentiles, output_address=mask_address)
obj.plot.cube('H1_6563A_b', Halpha_band, wcs=wcs, masks_file=mask_address)

# shoc579.fit.spatial_mask(spatial_mask_file, fit_conf=obs_cfg, line_detection=True, output_log=output_lines_log_file)


# obj.spatial_masker()

#             flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
#             flux_err = np.sqrt(cube[:, idx_j, idx_i].var.data) * norm_flux

# for i, obj in enumerate(objList):
#
#     # Data location
#     cube_address = fitsFolder/fileList[i]
#     objFolder = resultsFolder/obj
#     voxelFolder = resultsFolder/obj/'voxel_data'
#     db_addresss = objFolder/f'{obj}_database.fits'
#     maskFits_address = objFolder/f'{obj}_masks.fits'
#
#     # Output data:
#     fitsLog_addresss = objFolder/f'{obj}_linesLog.fits'
#     hdul_lineslog = fits.HDUList()
#
#     # Load data
#     wave, cube, header = import_muse_fits(cube_address)
#
#     # Loop throught the line regions
#
#     for idx_region in [0]:
#
#         # Voxel mask
#         start = time.time()
#         n_lines = 0
#         region_label = f'mask_{idx_region}'
#         region_mask = fits.getdata(maskFits_address, region_label, ver=1)
#         region_mask = region_mask.astype(bool)
#         idcs_voxels = np.argwhere(region_mask)
#
#         # Lines mask
#         mask_address = dataFolder/f'{obj}_region{idx_region}_mask.txt'
#         mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)
#         user_conf = obsData[f'region{idx_region}_line_fitting']
#
#         print(f'\n- Treating {region_label} with {idcs_voxels.shape[0]} pixels')
#
#         # Loop through voxels
#         n_voxels = idcs_voxels.shape[0]
#
#         # for idx_voxel in progressbar(np.arange(n_voxels), redirect_stdout=True):
#         for idx_voxel in np.arange(n_voxels):
#
#             idx_j, idx_i = idcs_voxels[idx_voxel]
#             voxel_dict = {}
#
#             grid_address_i = voxelFolder/f'{idx_j}-{idx_i}_LineGrid_{obj}.png'
#             pdfTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable'
#             txtTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable.txt'
#
#             flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
#             flux_err = np.sqrt(cube[:, idx_j, idx_i].var.data) * norm_flux
#
#             voxel = lime.Spectrum(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], norm_flux=norm_flux)
#
#             if verbose:
#                 voxel.plot_spectrum(spec_label=f'{obj} voxel {idx_j}-{idx_i}', log_scale=True)
#
#             # Security check for pixels with nan values:
#             idcs_nan = np.isnan(voxel.flux)
#             flux_interpolated = None
#
#             # Nans in flux
#             idcs_nan = np.isnan(voxel.flux)
#             if idcs_nan.any():
#                 Interpolation = interp1d(voxel.wave[~idcs_nan], voxel.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
#                 voxel.flux = Interpolation(voxel.wave)
#                 if verbose:
#                     voxel.plot_spectrum(comp_array=Interpolation(voxel.wave))
#
#             # Nans in errflux
#             idcs_nan = np.isnan(voxel.err_flux)
#             if idcs_nan.any():
#                 Interpolation = interp1d(voxel.wave[~idcs_nan], voxel.err_flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
#                 voxel.err_flux = Interpolation(voxel.wave)
#
#             peaks_table, matched_masks_DF = voxel.match_line_mask(mask_df, noise_region, detect_threshold=thres_array[idx_region],
#                                                                   width_mode='fixed')
#
#             # Index of detected lines
#             idcsObsLines = matched_masks_DF.index
#
#             if verbose:
#                 voxel.plot_spectrum(peaks_table=peaks_table, match_log=matched_masks_DF,
#                                     spec_label=f'{obj} voxel {idx_j}-{idx_i}', log_scale=True)
#
#             # Fit and check the regions
#             obsLines = matched_masks_DF.loc[idcsObsLines].index.values
#             for j, lineLabel in enumerate(obsLines):
#                 # if lineLabel == 'He1_5876A':
#                     wave_regions = matched_masks_DF.loc[lineLabel, 'w1':'w6'].values
#                     voxel.fit_from_wavelengths(lineLabel, wave_regions, user_cfg=user_conf)
#                     voxel.display_results()
#                 # try:
#                 #     voxel.fit_from_wavelengths(lineLabel, wave_regions, user_cfg=user_conf)
#                 #
#                 # except ValueError as e:
#                 #     err_value = 'NAN values' if 'NaN' in str(e) else 'valueError'
#                 #     err_label = f'ER_{lineLabel[lineLabel.find("_")+1:]}'
#                 #     voxel_dict[err_label] = err_value
#                 #     dict_errs[f'{idx_j}-{idx_i}_{lineLabel}'] = e
#                 #     print(f'--- Line measuring failure at {lineLabel} ({err_value}), {idx_j}-{idx_i}')
#
#             # Spectrum data
#             n_lines += len(voxel.log.index)
#             voxel_dict['N_Lines'] = len(voxel.log.index)
#             voxel_dict['N_nan'] = idcs_nan.sum()
#
#             # Converting linesLog to fits
#             linesHDU = lime.io.log_to_HDU(voxel.log, ext_name=f'{idx_j}-{idx_i}_linelog', header_dict=voxel_dict)
#
#             # Save spectrum data:
#             if linesHDU is not None:
#                 hdul_lineslog.append(linesHDU)
#
#         # Store the drive
#         # hdul_lineslog.writeto(fitsLog_addresss, overwrite=True, output_verify='fix')
#         end = time.time()
#         print(f'- Execution time {(end - start):.2f} seconds, {(end - start) / 60:.1f} min, for {n_lines} lines, voxels {n_voxels}')
#         # print(
#         #     f'- Execution time {(end - start):.2f} seconds, {(end - start) / 60:.1f} min, for {n_lines} lines, voxels {n_voxels},'
#         #     f' errors {len(dict_errs.keys())}')
#
# # Show summary
# print('-- Error summary')
# for voxel_fail, error in dict_errs.items():
#     print(voxel_fail)
