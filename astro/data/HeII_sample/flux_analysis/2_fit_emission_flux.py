import pathlib
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

conf_file_address = '../sampleHeII.ini'
obsData = sr.loadConfData(conf_file_address)

fits_folder = pathlib.Path(obsData['data_location']['fits_folder'])
treatment_folder = pathlib.Path(obsData['data_location']['treatment_folder'])
results_folder = pathlib.Path(obsData['data_location']['results_folder'])

catalogueDataframe = treatment_folder/f'AVO_catalogue_dataframe.txt'
normFlux = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']

logDF = sr.lineslogFile_to_DF(catalogueDataframe)
idcs_obj = logDF.valid == True
objList = logDF.loc[idcs_obj].index.values


mask_global = treatment_folder/'AVO_global_mask.txt'
maskDF = pd.read_csv(mask_global, delim_whitespace=True, header=0, index_col=0)


verbose = False

linesForced = ['He2_4686A', 'O3_4959A', 'O3_5007A', 'H1_4861A', 'H1_6563A', 'N2_6584A', 'S2_6716A', 'S2_6731A']

columns = ['Eqw_Halpha', 'Eqw_Halpha_err']
properties_df = pd.DataFrame(index=objList, columns=columns)
properties_df_addresss = treatment_folder/'AVO_properties_log.txt'
for line in linesForced:
    properties_df[f'{line}'] = np.nan
    properties_df[f'{line}_err'] = np.nan

print(f'Treating {objList.size}')
for i, obj in enumerate(objList):

    # if obj == '290-51941-177':

        # Inputs location
        print(f'\n- Treating {obj}')
        fits_address = fits_folder/f'{obj}.fits'
        objFolder = results_folder/'line_measurements'/obj

        # Outputs location
        local_mask = objFolder/f'mask_{obj}.txt'
        local_lineslog = objFolder/f'lineslog_{obj}.txt'
        pdf_lineslog = objFolder/f'tablelog_{obj}'

        # Make folder if not available
        objFolder.mkdir(parents=True, exist_ok=True)

        # Read the spectrum
        wave, data, hdrs = sr.import_fits_data(fits_address, instrument='SDSS')
        flux = data['flux'] * normFlux
        z_i = hdrs[1]["z"][0]

        lm = sr.LineMesurer(wave, flux, redshift=z_i, normFlux=normFlux)

        # ------------------------------ Check available lines
        idcs_nan = np.isnan(lm.flux)
        flux_interpolated = None

        if idcs_nan.any():
            Interpolation = interp1d(lm.wave[~idcs_nan], lm.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
            flux_interpolated = Interpolation(lm.wave)
            lm.flux = flux_interpolated
            norm_spec = lm.continuum_remover(noise_region)
        else:
            try:
                norm_spec = lm.continuum_remover(noise_region)
            except:
                norm_spec = None

        if norm_spec is not None:

            # Identify the emission lines
            obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=noise_region, intLineThreshold=1)
            local_maskDF = lm.match_lines(obsLinesTable, maskDF)

            # Add the following lines if they were missed:
            for lineLabel in linesForced:
                if local_maskDF.loc[lineLabel, 'observation'] == 'not detected':
                    local_maskDF.loc[lineLabel, 'observation'] = 'detected'

            idcsObsLines = (local_maskDF.observation == 'detected')
            if verbose:
                lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=local_maskDF, specLabel=f'{obj}')
                lm.plot_detected_lines(local_maskDF[idcsObsLines], ncols=8)


            # ------------------------------ Measure the fluxes
            lm = sr.LineMesurer(wave, flux, redshift=z_i, normFlux=normFlux)

            # Fit and check the regions
            obsLines = local_maskDF.loc[idcsObsLines].index.values
            for j, lineLabel in enumerate(obsLines):
                print(f'- {lineLabel}')
                wave_regions = local_maskDF.loc[lineLabel, 'w1':'w6'].values
                try:
                    lm.fit_from_wavelengths(lineLabel, wave_regions)
                    if lineLabel in ('He2_4686A', 'H1_6563A'):
                        # if verbose and False:
                        lm.plot_fit_components(lmfit_output=lm.fit_output, logScale=True, output_address=objFolder/f'{obj}_{lineLabel}_plot.png')
                except:
                    print(f'- Failure at: {lineLabel}')

            # Check Extinction
            lm.plot_line_grid(lm.linesDF, output_address=objFolder / f'{obj}_grid_plot.png')
            lm.save_lineslog(local_maskDF.loc[local_maskDF['observation'] == 'detected'], local_mask)
            lm.save_lineslog(lm.linesDF, local_lineslog)
            lm.table_fluxes(lm.linesDF, pdf_lineslog)

        for param in ('eqw', 'eqw_err'):
            if (param in lm.linesDF) and ('H1_6563A' in lm.linesDF.index):
                if param == 'eqw':
                    properties_df.loc[obj, 'Eqw_Halpha'] = lm.linesDF.loc['H1_6563A', param]
                else:
                    properties_df.loc[obj, 'Eqw_Halpha_err'] = lm.linesDF.loc['H1_6563A', param]

        for lineLabel in linesForced:
            if lineLabel in lm.linesDF.index:
                properties_df.loc[obj, f'{lineLabel}'] = lm.linesDF.loc[lineLabel, 'intg_flux']
                properties_df.loc[obj, f'{lineLabel}_err'] = lm.linesDF.loc[lineLabel, 'intg_err']

sr.save_lineslog(properties_df, properties_df_addresss)
            # if verbose:

                    # Save spectrum data:
                    # for key, value in voxel_dict.items():
                    #     linesHDU.header[key] = value
                    # hdul_lineslog.append(linesHDU)
                    # cHbeta, rc_pyneb = red_corr_HalphaHbeta_ratio(lm.linesDF, 0.0)
                    # lm.save_lineslog(lm.linesDF, local_lineslog)
                    # lm.table_fluxes(lm.linesDF, pdfTableFile, txtTableFile, rc_pyneb)

            # # Store the drive
            # hdul_lineslog.writeto(fitsLog_addresss, overwrite=True, output_verify='fix')
            # end = time.time()
            #
            # # Show summary
            # for voxel_fail, error in dict_errs.items():
            #     print(voxel_fail)
            # print(f'- Execution time {end - start:.3f}s, for {n_lines} lines, errors {len(dict_errs.keys())}')






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
#             local_maskDF = lm.match_lines(obsLinesTable, mask_df)
#             # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=local_maskDF)
#             idcsObsLines = (local_maskDF.observation == 'detected')
#             lm.plot_detected_lines(local_maskDF[idcsObsLines], ncols=8)
#
#             # Save the local mask dataframe
#             with open(local_mask, 'wb') as output_db:
#                 string_DF = local_maskDF.loc[idcsObsLines].to_string()
#                 output_db.write(string_DF.encode('UTF-8'))
#
#             # Reset and measure the lines
#             lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)
#
#             # Fit and check the regions
#             obsLines = local_maskDF.loc[idcsObsLines].index.values
#             for j, lineLabel in enumerate(obsLines):
#                 print(f'--- {lineLabel}:')
#                 wave_regions = local_maskDF.loc[lineLabel, 'w1':'w6'].values
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
#             lm.plot_detected_lines(local_maskDF[idcsObsLines], ncols=8, output_address=grid_address_i)
#
# print(dict_errs)
# try:
#     fits.update(fits_address, data=hdu.data, header=hdu.header, extname=extname)
# except:
#     fits.append(fits_address, data=hdu.data, header=hdu.header, extname=extname)