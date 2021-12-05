import numpy as np
import pandas as pd
import time
import lime
import progressbar
from scipy.interpolate import interp1d
from pathlib import Path
from astropy.io import fits


from src.specsiser.physical_model.extinction_model import ExtinctionModel
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits
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

    # Extinction model
    red_model = ExtinctionModel(Rv=obsData['Extinction']['R_v'], red_curve=obsData['Extinction']['red_law'])

    # Loop throught the line regions
    start = time.time()

    for idx_region in [0, 1, 2, 3, 4, 5]:

        # Voxel mask
        region_label = f'region_{idx_region}'
        region_mask = fits.getdata(maskFits_address, region_label, ver=1)
        region_mask = region_mask.astype(bool)
        idcs_voxels = np.argwhere(region_mask)

        # Lines mask
        mask_address = dataFolder/f'{obj}_region{idx_region}_mask.txt'
        mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)
        user_conf = obsData[f'region{idx_region}_line_fitting']

        print(f'\n- Treating {region_label} with {idcs_voxels.shape[0]} pixels')

        # Loop through voxels
        n_lines = 0
        n_voxels = idcs_voxels.shape[0]
        widgets = ['Voxels treated: ', progressbar.AnimatedMarker()]
        # bar = progressbar.ProgressBar(widgets=widgets).start()
        bar = progressbar.ProgressBar(maxval=n_voxels, widgets=widgets).start()
        for idx_voxel in np.arange(n_voxels):

            idx_j, idx_i = idcs_voxels[idx_voxel]
            voxel_dict = {}

            local_mask = voxelFolder/f'{idx_j}-{idx_i}_mask_{obj}.txt'
            local_lineslog = voxelFolder/f'{idx_j}-{idx_i}_lineslog_{obj}.txt'
            grid_address_i = voxelFolder/f'{idx_j}-{idx_i}_LineGrid_{obj}.png'
            pdfTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable'
            txtTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable.txt'

            flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
            flux_err = np.sqrt(cube[:, idx_j, idx_i].var.data) * norm_flux

            voxel = lime.Spectrum(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], norm_flux=norm_flux)

            if verbose:
                voxel.plot_spectrum(specLabel=f'{obj} voxel {idx_j}-{idx_i}', log_scale=True)

            # Security check for pixels with nan values:
            idcs_nan = np.isnan(voxel.flux)
            flux_interpolated = None

            if idcs_nan.any():
                if region_mask[idx_j, idx_i]:
                    Interpolation = interp1d(voxel.wave[~idcs_nan], voxel.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
                    if verbose:
                        voxel.plot_spectrum(continuumFlux=Interpolation(voxel.wave))
                    flux_interpolated = Interpolation(voxel.wave)
                    voxel.flux = flux_interpolated
                    norm_spec = lime.continuum_remover(voxel.wave_rest, voxel.flux, noise_region)
            else:
                norm_spec = lime.continuum_remover(voxel.wave_rest, voxel.flux, noise_region)

            # Identify the emission lines
            obsLinesTable = lime.line_finder(voxel.wave, norm_spec, noiseWaveLim=noise_region, intLineThreshold=3)
            maskLinesDF = lime.match_lines(voxel.wave_rest, voxel.flux, obsLinesTable, mask_df)
            idcsObsLines = (maskLinesDF.observation == 'detected')

            if verbose:
                voxel.plot_spectrum(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF, specLabel=f'{obj} voxel {idx_j}-{idx_i}')

            # Reset and measure the lines
            voxel = lime.Spectrum(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], norm_flux=norm_flux)

            # For central pixels with nan entries
            idcs_nan = np.isnan(voxel.flux)
            if idcs_nan.any():
                if region_mask[idx_j, idx_i]:
                    Interpolation = interp1d(voxel.wave[~idcs_nan], voxel.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
                    Interpolation_err = interp1d(voxel.wave[~idcs_nan], voxel.errFlux[~idcs_nan], kind='slinear', fill_value="extrapolate")
                    voxel.flux = Interpolation(voxel.wave)
                    voxel.errFlux = Interpolation_err(voxel.wave)

            # Fit and check the regions
            obsLines = maskLinesDF.loc[idcsObsLines].index.values
            for j, lineLabel in enumerate(obsLines):
                wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
                try:
                    voxel.fit_from_wavelengths(lineLabel, wave_regions, user_cfg=user_conf)
                    if verbose:
                        voxel.plot_fit_components(voxel.fit_output, log_scale=True)

                except ValueError as e:
                    err_value = 'NAN values' if 'NaN' in str(e) else 'valueError'
                    err_label = f'ER_{lineLabel[lineLabel.find("_")+1:]}'
                    voxel_dict[err_label] = err_value
                    dict_errs[f'{idx_j}-{idx_i}_{lineLabel}'] = e
                    print(f'--- Line measuring failure at {lineLabel} ({err_value}), {idx_j}-{idx_i}')

            # Check Extinction
            if verbose:
                voxel.plot_line_grid(voxel.linesDF)
                cHbeta, cHbeta_err = red_model.cHbeta_from_log(voxel.linesDF, plot_address=True)

            # Spectrum data
            n_lines += len(voxel.linesDF.index)
            voxel_dict['N_Lines'] = len(voxel.linesDF.index)
            voxel_dict['N_nan'] = idcs_nan.sum()

            # Converting linesLog to fits
            linesHDU = lime.io.lineslog_to_HDU(voxel.linesDF, ext_name=f'{idx_j}-{idx_i}_linelog', header_dict=voxel_dict)

            # Save spectrum data:
            hdul_lineslog.append(linesHDU)
            bar.update(idx_voxel)

        # Store the drive
        hdul_lineslog.writeto(fitsLog_addresss, overwrite=True, output_verify='fix')
        end = time.time()

# Show summary
for voxel_fail, error in dict_errs.items():
    print(voxel_fail)
print(f'- Execution time {(end - start)/60:.3f} min, for {n_lines} lines, errors {len(dict_errs.keys())}')






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