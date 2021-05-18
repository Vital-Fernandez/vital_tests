import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, red_corr_HalphaHbeta_ratio, store_frame_to_fits
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
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
        flux6312_levels = np.nanpercentile(flux6312_image, pertil_array)

        flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
        flux6563_levels = np.nanpercentile(flux6563_image, pertil_array)

        Halpha_min_level = flux6563_levels[3]
        sulfur_levels = np.array([0, 1, 2])

        for idx_sulfur in sulfur_levels:

            # Search within that limit
            if idx_sulfur == 0:
                maFlux_image = np.ma.masked_where((flux6312_image >= flux6312_levels[idx_sulfur]) &
                                                  (flux6563_image > Halpha_min_level),
                                                  flux6563_image)
            else:
                maFlux_image = np.ma.masked_where((flux6312_image >= flux6312_levels[idx_sulfur]) &
                                                  (flux6312_image < flux6312_levels[idx_sulfur-1]) &
                                                  (flux6563_image > Halpha_min_level),
                                                  flux6563_image)
                print(flux6312_levels[idx_sulfur], flux6312_levels[idx_sulfur-1])


            idcs_voxels = np.argwhere(maFlux_image.mask)
            print(f'({idcs_voxels.shape[0]} pixels)')

            # if verbose:
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
            ax.update({'title': r'{} galaxy, $H\alpha$ flux'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
            im = ax.imshow(maFlux_image, cmap=cm.gray,
                           norm=colors.SymLogNorm(linthresh=flux6563_levels[-2], vmin=flux6563_levels[-2], base=10))
            plt.show()


            mask_name = f'region_{idx_sulfur}'
            mask_hdu = fits.ImageHDU(name=mask_name, data=maFlux_image.mask.astype(int), ver=1)
            store_frame_to_fits(db_addresss, fits_hdu=mask_hdu, ext_name=mask_name)


        # --------------------  Clusters mask
        maFlux_image = np.ma.masked_where((flux6312_image < flux6312_levels[sulfur_levels[-1]]) &
                                          (flux6563_image >= Halpha_min_level),
                                          flux6563_image)

        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
        ax.update({'title': r'{} galaxy, $H\alpha$ flux'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
        im = ax.imshow(maFlux_image, cmap=cm.gray,
                       norm=colors.SymLogNorm(linthresh=flux6563_levels[-2], vmin=flux6563_levels[-2], base=10))
        plt.show()

        mask_name = f'region_3'
        mask_hdu = fits.ImageHDU(name=mask_name, data=maFlux_image.mask.astype(int), ver=1)
        store_frame_to_fits(db_addresss, fits_hdu=mask_hdu, ext_name=mask_name)


        # ------------------------ Intermediate galaxy
        maFlux_image = np.ma.masked_where((flux6563_image < Halpha_min_level) &
                                          (flux6563_image >= flux6563_levels[4]),
                                          flux6563_image)

        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
        ax.update({'title': r'{} galaxy, $H\alpha$ flux'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
        im = ax.imshow(maFlux_image, cmap=cm.gray,
                       norm=colors.SymLogNorm(linthresh=flux6563_levels[-2], vmin=flux6563_levels[-2], base=10))
        plt.show()

        mask_name = f'region_4'
        mask_hdu = fits.ImageHDU(name=mask_name, data=maFlux_image.mask.astype(int), ver=1)
        store_frame_to_fits(db_addresss, fits_hdu=mask_hdu, ext_name=mask_name)



        # ------------------------ Complete galaxy
        maFlux_image = np.ma.masked_where((flux6563_image < flux6563_levels[4]) &
                                          (flux6563_image >= flux6563_levels[5]),
                                          flux6563_image)

        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
        ax.update({'title': r'{} galaxy, $H\alpha$ flux'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
        im = ax.imshow(maFlux_image, cmap=cm.gray,
                       norm=colors.SymLogNorm(linthresh=flux6563_levels[-2], vmin=flux6563_levels[-2], base=10))
        plt.show()

        mask_name = f'region_5'
        mask_hdu = fits.ImageHDU(name=mask_name, data=maFlux_image.mask.astype(int), ver=1)
        store_frame_to_fits(db_addresss, fits_hdu=mask_hdu, ext_name=mask_name)


        # print(f'\nUsing line Halpha at percentile {pertil_array[hydrogen_bdry]} = {flux6563_levels[hydrogen_bdry]:.2f}'
        #       f' ({idcs_voxels.shape[0]} pixels)')

        # # Loop through voxels
        # n_lines = 0
        # start = time.time()
        # for idx_voxel, idx_pair in enumerate(idcs_voxels):
        #
        #     print(f'-- Treating voxel {idx_voxel} {idx_pair}')
        #     idx_j, idx_i = idx_pair
        #     voxel_dict = {}
        #
        #     local_mask = voxelFolder/f'{idx_j}-{idx_i}_mask_{obj}.txt'
        #     local_lineslog = voxelFolder/f'{idx_j}-{idx_i}_lineslog_{obj}.txt'
        #     grid_address_i = voxelFolder/f'{idx_j}-{idx_i}_LineGrid_{obj}.png'
        #     pdfTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable'
        #     txtTableFile = voxelFolder/f'{idx_j}-{idx_i}_linesTable.txt'
        #
        #     flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
        #     flux_err = cube[:, idx_j, idx_i].var.data * norm_flux
        #
        #     lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)
        #
        #     if verbose:
        #         lm.plot_spectrum_components(specLabel=f'{obj} voxel {idx_j}-{idx_i}', log_scale=True)
        #
        #     # Security check for pixels with nan values:
        #     idcs_nan = np.isnan(lm.flux)
        #
        #     if idcs_nan.any():
        #         Interpolation = interp1d(lm.wave[~idcs_nan], lm.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
        #         lm.flux = Interpolation(lm.wave)
        #         norm_spec = lm.continuum_remover(noise_region)
        #     else:
        #         norm_spec = lm.continuum_remover(noise_region)
        #
        #     # Identify the emission lines
        #     obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=noise_region, intLineThreshold=3)
        #     maskLinesDF = lm.match_lines(obsLinesTable, mask_df)
        #     idcsObsLines = (maskLinesDF.observation == 'detected')
        #
        #     if verbose:
        #         lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF, specLabel=f'{obj} voxel {idx_j}-{idx_i}')
        #         lm.plot_detected_lines(maskLinesDF[idcsObsLines], ncols=8)
        #
        #     # Reset and measure the lines
        #     lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)
        #
        #     # Fit and check the regions
        #     obsLines = maskLinesDF.loc[idcsObsLines].index.values
        #     for j, lineLabel in enumerate(obsLines):
        #         wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
        #         try:
        #             lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf={})
        #         except ValueError as e:
        #             err_value = 'NAN values' if 'NaN' in str(e) else 'valueError'
        #             err_label = f'ER_{lineLabel[lineLabel.find("_")+1:]}'
        #             voxel_dict[err_label] = err_value
        #             dict_errs[f'{idx_j}-{idx_i}_{lineLabel}'] = e
        #             print(f'--- Line measuring failure at {lineLabel} ({err_value})')
        #
        #     # Check Extinction
        #     if verbose:
        #         print('Show extinction calculation')
        #         lm.plot_line_grid(lm.linesDF)
        #         cHbeta, cHbeta_err = red_model.cHbeta_from_log(lm.linesDF, plot_address=True)
        #
        #     # Spectrum data
        #     n_lines += len(lm.linesDF.index)
        #     voxel_dict['N_Lines'] = len(lm.linesDF.index)
        #     voxel_dict['N_nan'] = idcs_nan.sum()
        #
        #     # Converting linesLog to fits
        #     linesSA = lm.linesDF.to_records(index=True, column_dtypes=default_linelog_types, index_dtypes='<U50')
        #     linesCol = fits.ColDefs(linesSA)
        #     linesHDU = fits.BinTableHDU.from_columns(linesCol, name=f'{idx_j}-{idx_i}_linelog')
        #
        #     # Save spectrum data:
        #     for key, value in voxel_dict.items():
        #         linesHDU.header[key] = value
        #     hdul_lineslog.append(linesHDU)
        #     # cHbeta, rc_pyneb = red_corr_HalphaHbeta_ratio(lm.linesDF, 0.0)
        #     # lm.save_lineslog(lm.linesDF, local_lineslog)
        #     # lm.table_fluxes(lm.linesDF, pdfTableFile, txtTableFile, rc_pyneb)
        #
        # # Store the drive
        # # hdul_lineslog.writeto(fitsLog_addresss, overwrite=True, output_verify='fix')
        # end = time.time()
        #
        # # Show summary
        # for voxel_fail, error in dict_errs.items():
        #     print(voxel_fail)
        # print(f'- Execution time {end - start:.3f}s, for {n_lines} lines, errors {len(dict_errs.keys())}')