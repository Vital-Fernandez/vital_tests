import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, red_corr_HalphaHbeta_ratio, default_linelog_types
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams
from astropy.wcs import WCS
import time

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

dict_errs = {}
dict_nan_values = {}

ref_flux_line = 'S3_6312A'
int_flux_boundary = -3

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        cube_address = fitsFolder/fileList[i]
        objFolder = resultsFolder/obj
        voxelFolder = resultsFolder/obj/'voxel_data'
        db_addresss = objFolder / f'{obj}_database.fits'
        mask_address = dataFolder/obj/f'{obj}_mask.txt'

        # Output data:
        fitsLog_addresss = objFolder / f'{obj}_linesLog.fits'
        hdul_lineslog = fits.HDUList()

        # Load data
        wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')
        wave_rest = wave / (1 + z_objs[i])
        mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)

        # Declare voxels to analyse
        flux6312_image = fits.getdata(db_addresss, f'{ref_flux_line}_flux', ver=1)
        flux6312_contours = fits.getdata(db_addresss, f'{ref_flux_line}_contour', ver=1)
        flux6312_levels = np.percentile(flux6312_contours[flux6312_contours > 0], pertil_array)
        flux6312_boundary = flux6312_levels[int_flux_boundary]

        flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
        flux6563_contours = fits.getdata(db_addresss, f'H1_6563A_contour', ver=1)
        flux6563_levels = np.percentile(flux6563_contours[flux6563_contours > 0], pertil_array)
        flux6563_boundary = flux6563_levels[int_flux_boundary]

        idcs_bolean = (flux6312_contours >= flux6312_boundary) & (flux6563_contours >= flux6563_boundary)
        idcs_voxels = np.argwhere(idcs_bolean)
        print(f'\nUsing line {ref_flux_line} at percentile {pertil_array[int_flux_boundary]} = {flux6312_boundary:.2f} ({idcs_voxels.shape[0]} pixels)')

        # Loop through voxels
        n_lines = 0
        start = time.time()
        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            print(f'-- Treating voxel {idx_voxel} {idx_pair}')
            idx_j, idx_i = idx_pair
            voxel_dict = {}

            local_mask = voxelFolder/f'{idx_j}-{idx_i}_mask_{obj}.txt'
            local_lineslog = voxelFolder/f'{idx_j}-{idx_i}_lineslog_{obj}.txt'
            grid_address_i = voxelFolder/f'{idx_j}-{idx_i}_LineGrid_{obj}.png'
            pdfTableFile = voxelFolder / f'{idx_j}-{idx_i}_linesTable'
            txtTableFile = voxelFolder / f'{idx_j}-{idx_i}_linesTable.txt'

            flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
            flux_err = cube[:, idx_j, idx_i].var.data * norm_flux

            lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)
            # lm.plot_spectrum_components()

            # Security check for pixels with nan values:
            idcs_nan = np.isnan(lm.flux)

            if idcs_nan.any():
                Interpolation = interp1d(lm.wave[~idcs_nan], lm.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
                lm.flux = Interpolation(lm.wave)
                norm_spec = lm.continuum_remover(noise_region)
            else:
                norm_spec = lm.continuum_remover(noise_region)

            # Identify the emission lines
            obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=noise_region, intLineThreshold=3)
            maskLinesDF = lm.match_lines(obsLinesTable, mask_df)
            # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF)
            idcsObsLines = (maskLinesDF.observation == 'detected')
            # lm.plot_detected_lines(maskLinesDF[idcsObsLines], ncols=8)

            # Reset and measure the lines
            lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)

            # Fit and check the regions
            obsLines = maskLinesDF.loc[idcsObsLines].index.values
            for j, lineLabel in enumerate(obsLines):
                wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
                try:
                    lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf={})
                except ValueError as e:
                    err_value = 'NAN values' if 'NaN' in str(e) else 'valueError'
                    err_label = f'ER_{lineLabel[lineLabel.find("_")+1:]}'
                    voxel_dict[err_label] = err_value
                    print(f'--- Line measuring failure at {lineLabel} ({err_value})')

            # Spectrum data
            n_lines += len(lm.linesDF.index)
            voxel_dict['N_Lines'] = len(lm.linesDF.index)
            voxel_dict['N_nan'] = idcs_nan.sum()

            # Converting linesLog to fits
            linesSA = lm.linesDF.to_records(index=True, column_dtypes=default_linelog_types, index_dtypes='<U50')
            linesCol = fits.ColDefs(linesSA)
            linesHDU = fits.BinTableHDU.from_columns(linesCol, name=f'{idx_j}-{idx_i}_linelog')

            # Save spectrum data:
            for key, value in voxel_dict.items():
                linesHDU.header[key] = value
            hdul_lineslog.append(linesHDU)
            # cHbeta, rc_pyneb = red_corr_HalphaHbeta_ratio(lm.linesDF, 0.0)
            # lm.save_lineslog(lm.linesDF, local_lineslog)
            # lm.table_fluxes(lm.linesDF, pdfTableFile, txtTableFile, rc_pyneb)

        # Store the drive
        hdul_lineslog.writeto(fitsLog_addresss, overwrite=True, output_verify='fix')
        end = time.time()
        print(f'- Execution time {end - start:.3f}s, for {n_lines} lines')




#             print(f'- Printing results tables')
#             lm.save_lineslog(lm.linesDF, local_lineslog)
#             lm.table_fluxes(lm.linesDF, pdfTableFile, txtTableFile, rc_pyneb)
#             lm.plot_detected_lines(maskLinesDF[idcsObsLines], ncols=8, output_address=grid_address_i)
#
# print(dict_errs)

#       IT IS NOT A BAD CODE
#         # Loop through voxels
#         for idx_voxel, idx_pair in enumerate(idcs_voxels):
#
#             print(f'-- Treating voxel {idx_voxel} {idx_pair}')
#             idx_j, idx_i = idx_pair
#             voxel_dict = {}
#
#             idx_database = (obj_db.y_voxel == idx_j) & (obj_db.x_voxel == idx_i)
#             local_mask = voxelFolder/f'{idx_j}-{idx_i}_mask_{obj}.txt'
#             local_lineslog = voxelFolder/f'{idx_j}-{idx_i}_lineslog_{obj}.txt'
#             grid_address_i = voxelFolder/f'{idx_j}-{idx_i}_LineGrid_{obj}.png'
#             pdfTableFile = voxelFolder / f'{idx_j}-{idx_i}_linesTable'
#             txtTableFile = voxelFolder / f'{idx_j}-{idx_i}_linesTable.txt'
#
#             flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
#             flux_err = cube[:, idx_j, idx_i].var.data * norm_flux
#
#             lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)
#             lm.plot_spectrum_components()
#
#             # Normalize
#             norm_spec = lm.continuum_remover(noise_region)
#
#             # Security check for pixels with nan values:
#             idcs_nan = np.isnan(lm.flux)
#
#             if idcs_nan.any():
#                 Interpolation = interp1d(lm.wave[~idcs_nan], lm.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
#                 lm.flux = Interpolation(lm.wave)
#                 norm_spec = lm.continuum_remover(noise_region)
#                 #obj_db.loc[idx_database, 'Bad_values'] = idcs_nan.sum()
#
#             # Identify the emission lines
#             obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=noise_region, intLineThreshold=3)
#             maskLinesDF = lm.match_lines(obsLinesTable, mask_df)
#             # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF)
#             idcsObsLines = (maskLinesDF.observation == 'detected')
#             lm.plot_detected_lines(maskLinesDF[idcsObsLines], ncols=8)
#
#             # Save the local mask dataframe
#             with open(local_mask, 'wb') as output_db:
#                 string_DF = maskLinesDF.loc[idcsObsLines].to_string()
#                 output_db.write(string_DF.encode('UTF-8'))
#
#             # Reset and measure the lines
#             lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)
#
#             # Fit and check the regions
#             obsLines = maskLinesDF.loc[idcsObsLines].index.values
#             for j, lineLabel in enumerate(obsLines):
#                 print(f'--- {lineLabel}:')
#                 wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
#                 try:
#                     lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf={})
#                     print(lm)
#                     lm.plot_fit_components(lm.fit_output)
#                 except:
#                     if lineLabel == 'H1_6563A':
#                         obj_db.loc[idx_database, 'Halpha_nan_pixel'] = True
#                     else:
#                         dict_errs[f'{lineLabel}_{idx_database}'] = 'Err'
#
#             # Number of lines measured
#             obj_db.loc[idx_database, 'n_emissions'] = len(lm.linesDF.index)
#
#             # Storing special fluxes in database
#             for label in columns_fluxes:
#                 if label in lm.linesDF.index:
#                     obj_db.loc[idx_database, label] = lm.linesDF.loc[label, 'intg_flux']
#                     obj_db.loc[idx_database, label + '_err'] = lm.linesDF.loc[label, 'intg_err']
#
#             # Compute reddening correction
#             cHbeta, rc_pyneb = red_corr_HalphaHbeta_ratio(lm.linesDF, 0.0)
#             obj_db.loc[idx_database, 'cHbeta'] = cHbeta
#
#             # Save the results
#             print(f'- Printing results tables')
#             lm.save_lineslog(lm.linesDF, local_lineslog)
#             lm.table_fluxes(lm.linesDF, pdfTableFile, txtTableFile, rc_pyneb)
#             lm.plot_detected_lines(maskLinesDF[idcsObsLines], ncols=8, output_address=grid_address_i)
#
# print(dict_errs)
# try:
#     fits.update(fits_address, data=hdu.data, header=hdu.header, extname=extname)
# except:
#     fits.append(fits_address, data=hdu.data, header=hdu.header, extname=extname)