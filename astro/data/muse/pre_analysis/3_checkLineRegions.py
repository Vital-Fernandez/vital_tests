import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from mpdaf.obj import Cube
import matplotlib.pyplot as plt
import astropy.units as u
from mpdaf.obj import deg2sexa
from astropy.wcs import WCS


# Declare data and files location
conf_file_address = '../muse_greenpeas.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)
objList = obsData['sample_data']['object_list']
fileList = obsData['sample_data']['file_list']
dataFolder = obsData['sample_data']['data_folder']
z_objs = obsData['sample_data']['z_array']
noise_region = obsData['sample_data']['noiseRegion_array']
mask_columns = ['wavelength', 'ion', 'w1', 'w2', 'w3', 'w4', 'w5', 'w6']

voxel_list = [(105, 130), (92, 102)]
window_box = 2

range_box = np.arange(-window_box, window_box + 1)


for idx_obj, obj in enumerate(objList):

    if idx_obj == 0:

        # Load the data
        print(f'\n- {obj}')
        fits_address_i = f'{dataFolder}/{fileList[idx_obj]}'
        wave, cube, header = sr.import_fits_data(fits_address_i, instrument='MUSE')
        wave_rest = wave / (1 + z_objs[idx_obj])
        center_voxel = voxel_list[idx_obj]

        # Loop through the voxels
        counter = 0

        for i in range_box:
            for j in range_box:
                idx_i, idx_j = center_voxel[0] + i, center_voxel[1] + j
                print(f'-- Voxel {counter}: {idx_i} {idx_j}')

                # Declare voxels
                voxel = cube[:, idx_i, idx_j]
                flux_voxel = voxel.data.data

                if not np.isnan(flux_voxel).any():

                    lm = sr.LineMesurer(wave_rest, flux_voxel)

                    # Identify the emission lines
                    norm_flux = lm.continuum_remover(noise_region)
                    obsLinesTable = lm.line_finder(norm_flux, noiseWaveLim=noise_region, intLineThreshold=3)
                    obsLinesDF = lm.match_lines(obsLinesTable, sr._linesDb)
                    lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=obsLinesDF)

                    # # Get matched lines
                    # idcsObsLines = (obsLinesDF.observation == 'detected')
                    # obsLines = obsLinesDF.loc[idcsObsLines].index.values
                    #
                    # # Fit and check the regions
                    # logs_name_i = fileList[idx_obj].replace(".fits", f"_{idx_i}-{idx_j}_lineLog.txt")
                    # lineslog_address_i = f'{dataFolder}/{logs_name_i}'
                    # lm = sr.LineMesurer(wave_rest, flux_voxel)
                    # for j, lineLabel in enumerate(obsLines):
                    #     print(f'-- {lineLabel}:')
                    #     wave_regions = obsLinesDF.loc[lineLabel, 'w1':'w6'].values
                    #     lm.fit_from_wavelengths(lineLabel, wave_regions)
                    # lm.save_lineslog(lm.linesDF, lineslog_address_i)
                    #
                    # # Plot checking model
                    # lm.linesLogAddress = lineslog_address_i
                    # lm.plot_detected_lines(lm.linesDF, ncols=8)

                else:
                    print('-- Warning NAN flux entries')

                counter += 1


        # # Load the data
        # print(f'\n- {obj}')
        # fits_address_i = f'{dataFolder}/{fileList[idx_obj]}'
        # mask_address_i = f'{dataFolder}/{fileList[idx_obj].replace(".fits", "_mask.txt")}'
        # lineslog_address_i = f'{dataFolder}/{fileList[idx_obj].replace(".fits", "_lineLog.txt")}'
        # linesGrid_address_i = f'{dataFolder}/{fileList[idx_obj].replace(".fits", "_lineGrid.png")}'
        # maskDF = pd.read_csv(mask_address_i, delim_whitespace=True, header=0, index_col=0)

        # # Declare spectrum
        # wave, cube, header = sr.import_fits_data(fits_address_i, instrument='MUSE')
        # voxel = cube[:, voxel_list[i][0], voxel_list[i][1]]
        # flux_voxel = voxel.data.data
        # wave_rest = wave / (1 + z_objs[i])
        #
        # # Loop through the lines
        # lm = sr.LineMesurer(wave_rest, flux_voxel, lineslog_address_i)
        # obsLines = maskDF.index.values
        # for j, lineLabel in enumerate(obsLines):
        #     # Fit each line regions data
        #     print(f'-- {lineLabel}:')
        #     wave_regions = maskDF.loc[lineLabel, 'w1':'w6'].values
        #     lm.fit_from_wavelengths(lineLabel, wave_regions)
        # lm.save_lineslog(lm.linesDF, lineslog_address_i)
        #
        # # Plot checking model
        # lm.plot_detected_lines(lm.linesDF, ncols=8)



        # # Plot the single lines:
        # idcs_unblended = ~lm.linesDF.index.str.contains('_b')
        # lm.plot_line_grid(lm.linesDF.loc[idcs_unblended], ncols=8, output_address=linesGrid_address_i)


        # lm.linesDF = pd.DataFrame(columns=sr._linesDb.columns.values, index=maskDF.index.values)
        # for column in maskDF:
        #     lm.linesDF[column] = maskDF[column]