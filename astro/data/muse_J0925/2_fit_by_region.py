import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import compute_line_flux_image, red_corr_HalphaHbeta_ratio

# Declare data and files location
obsData = sr.loadConfData('muse_J0925.ini', group_variables=False)
objList = np.array([obsData['sample_data']['object_list']])
fileList = np.array([obsData['sample_data']['file_list']])
dataFolder = Path(obsData['sample_data']['data_folder'])
resultsFolder = Path(obsData['sample_data']['results_folder'])
z_objs = np.array([obsData['sample_data']['z_array']])
pertil_array = obsData['sample_data']['percentil_array']
db_headers = obsData['sample_data']['database_header_list']
db_format = {'DEC_deg': '{: 0.8f}', 'RA_deg': '{: 0.8f}'}
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']
dict_errs = {}

lineAreas = {'H1_6563A_b': (6533.0, 6596.0),
             'S3_6312A': (6310.0, 6319.0),
             'O3_5007A': (4999.0, 5025.0),
             'O3_4363A': (4355.0, 4374.0)}

export_elements = ['intg_flux', 'intg_err', 'gauss_flux', 'gauss_err', 'amp', 'mu', 'sigma', 'amp_err', 'mu_err',
                   'sigma_err', 'v_r', 'v_r_err', 'sigma_vel', 'sigma_err_vel']

for i, obj in enumerate(objList):

        # Data location
        objFolder = resultsFolder
        db_addresss = resultsFolder/f'{obj}_database.txt'
        voxelFolder = resultsFolder/'voxel_data'
        cube_address = dataFolder/fileList[i]
        mask_address = objFolder/f'{obj}_mask.txt'

        # Load data
        wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')
        wave_rest = wave / (1 + z_objs[i])
        obj_db = pd.read_csv(db_addresss, delim_whitespace=True, header=0, index_col=0)
        mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)

        # New database columns
        obj_db['Bad_values'] = 0.0
        obj_db['Halpha_nan_pixel'] = False
        obj_db['n_emissions'] = 0.0
        obj_db['cHbeta'] = 0.0
        obj_db['cHbeta_err'] = 0.0

        for lineComp in obsData['default_line_fitting']['H1_6563A_b'].split('-'):
            for param in export_elements:
                column_name = f'{lineComp}-{param}'
                obj_db[column_name] = np.nan
                print(column_name)

        for label in lineAreas:
            obj_db[label] = 'none'
            obj_db[label+'_err'] = 'none'

        # For testing
        obj_db['myThings'] = np.nan
        obj_db['peak_max'] = np.nan

        print(f'\n- {obj}: Cube dimensions {cube.shape}')

        # Get line region data
        lineFlux_dict, levelFlux_dict, levelText_dict = compute_line_flux_image(lineAreas,
                                                                                cube,
                                                                                z_objs[i],
                                                                                percent_array=pertil_array)

        # Declare voxels to analyse
        fluxImage_5007 = lineFlux_dict['H1_6563A_b']
        fluxLevels_5007 = levelFlux_dict['H1_6563A_b']

        image_size = fluxImage_5007.shape
        y, x = np.ogrid[0:image_size[0], 0:image_size[1]]

        int_level = fluxLevels_5007[-4]#[-4]
        mask_flux = fluxImage_5007 >= int_level
        idcs_voxels = np.argwhere(mask_flux)
        n_voxels = idcs_voxels.size / 2
        np.max(fluxImage_5007[mask_flux])

        # Loop through voxels
        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            print(f'-- Treating voxel {idx_voxel} {idx_pair} of {n_voxels}')
            idx_j, idx_i = idx_pair

            idx_database = (obj_db.y_voxel == idx_j) & (obj_db.x_voxel == idx_i)

            local_mask = voxelFolder/f'{idx_j}-{idx_i}_mask_{obj}.txt'
            local_lineslog = voxelFolder/f'{idx_j}-{idx_i}_lineslog_{obj}.txt'
            grid_address_i = voxelFolder/f'{idx_j}-{idx_i}_LineGrid_{obj}.png'
            fitcomponents_Halpha = voxelFolder/f'{idx_j}-{idx_i}_Halpha_{obj}.png'
            voxelSpectrum = voxelFolder / f'{idx_j}-{idx_i}_voxel_spec_{obj}.png'
            pdfTableFile = voxelFolder / f'{idx_j}-{idx_i}_linesTable'
            txtTableFile = voxelFolder / f'{idx_j}-{idx_i}_linesTable.txt'

            flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
            flux_err = cube[:, idx_j, idx_i].var.data * norm_flux

            lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)
            lm.plot_spectrum_components(output_address=voxelSpectrum)

            # Normalize
            try:
                norm_spec = lm.continuum_remover(noise_region)
                security_check = True
            except:
                security_check = False

            if security_check:
                # Security check for pixels with nan values:
                idcs_nan = np.isnan(lm.flux)

                if idcs_nan.any():
                    Interpolation = interp1d(lm.wave[~idcs_nan], lm.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
                    lm.flux = Interpolation(lm.wave)
                    norm_spec = lm.continuum_remover(noise_region)
                    obj_db.loc[idx_database, 'Bad_values'] = idcs_nan.sum()

                # Identify the emission lines
                obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=noise_region, intLineThreshold=3)
                maskLinesDF = lm.match_lines(obsLinesTable, mask_df, find_line_borders=False)
                # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF)
                idcsObsLines = (maskLinesDF.observation == 'detected')
                # lm.plot_detected_lines(maskLinesDF[idcsObsLines], ncols=8)

                # Save the local mask dataframe
                with open(local_mask, 'wb') as output_db:
                    string_DF = maskLinesDF.loc[idcsObsLines].to_string()
                    output_db.write(string_DF.encode('UTF-8'))

                # Reset and measure the lines
                lm = sr.LineMesurer(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], normFlux=norm_flux)

                # Fit and check the regions
                obsLines = maskLinesDF.loc[idcsObsLines].index.values
                for j, lineLabel in enumerate(obsLines):
                    print(f'--- {lineLabel}:')
                    wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
                    try:
                        lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=obsData['default_line_fitting'])
                        # lm.print_results(show_fit_report=True, show_plot=True)

                        if lineLabel == 'H1_6563A_b':
                            lm.plot_fit_components(lm.fit_output, output_address=fitcomponents_Halpha)
                    except:
                        if lineLabel == 'H1_6563A_b':
                            obj_db.loc[idx_database, 'Halpha_nan_pixel'] = True
                        else:
                            dict_errs[f'{lineLabel}_{idx_database}'] = 'Err'

                if 'H1_6563A' in lm.linesDF.index:

                    # Number of lines measured
                    obj_db.loc[idx_database, 'n_emissions'] = len(lm.linesDF.index)

                    # Storing special fluxes in database
                    for lineComp in obsData['default_line_fitting']['H1_6563A_b'].split('-'):
                        for param in export_elements:
                            if lineComp in lm.linesDF.index:
                                column_name = f'{lineComp}-{param}'
                                obj_db.loc[idx_database, column_name] = lm.linesDF.loc[lineComp, param]

                    # Compute reddening correction
                    cHbeta, rc_pyneb = red_corr_HalphaHbeta_ratio(lm.linesDF, 0.0)
                    obj_db.loc[idx_database, 'cHbeta'] = cHbeta

                    # Save the results
                    print(f'- Printing results tables')
                    idcs_obsLines = ~lm.linesDF.index.str.contains('_b')
                    lm.save_lineslog(lm.linesDF, local_lineslog)
                    # lm.table_fluxes(lm.linesDF[idcs_obsLines], pdfTableFile, txtTableFile, rc_pyneb)
                    # lm.plot_detected_lines(maskLinesDF[idcsObsLines], ncols=8, output_address=grid_address_i)

        # Save the object database
        with open(db_addresss, 'wb') as output_db:
            string_DF = obj_db.to_string()
            output_db.write(string_DF.encode('UTF-8'))

print(dict_errs)
