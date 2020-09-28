import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from mpdaf.obj import Cube
import matplotlib.pyplot as plt
import astropy.units as u
from mpdaf.obj import deg2sexa
from astropy.wcs import WCS

def lineMeanSampleProperties(label, dictDF):
    w1, w2, w3, w4, w5, w6 = [], [], [], [], [], []
    for id_voxel, df_voxel in dictDF.items():
        if label in df_voxel.index:
            w1.append(df_voxel.loc[label, 'w1'])
            w2.append(df_voxel.loc[label, 'w2'])
            w3.append(df_voxel.loc[label, 'w3'])
            w4.append(df_voxel.loc[label, 'w4'])
            w5.append(df_voxel.loc[label, 'w5'])
            w6.append(df_voxel.loc[label, 'w6'])
    w1, w2, w3, w4, w5, w6 = np.array(w1), np.array(w2), np.array(w3), np.array(w4), np.array(w5), np.array(w6)
    print(f'-- {label}: {w1.mean():.2f} {w2.mean():.2f} {w3.mean():.2f} {w4.mean():.2f} {w5.mean():.2f} {w6.mean():.2f}')
    return

# Declare data and files location
conf_file_address = '../muse_greenpeas.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)
objList = obsData['sample_data']['object_list']
fileList = obsData['sample_data']['file_list']
dataFolder = obsData['sample_data']['data_folder']
z_objs = obsData['sample_data']['z_array']
noise_region = obsData['sample_data']['noiseRegion_array']
mask_columns = ['wavelength', 'ion', 'w1', 'w2', 'w3', 'w4', 'w5', 'w6']
voxel_list = [(170, 170), (92, 102)]
window_box = 2
range_box = np.arange(-window_box, window_box + 1)

df_dict = {}
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

                    logs_name_i = fileList[idx_obj].replace(".fits", f"_{idx_i}-{idx_j}_lineLog.txt")
                    lineslog_address_i = f'{dataFolder}/{logs_name_i}'
                    df_dict[f'{idx_i}-{idx_j}'] = pd.read_csv(lineslog_address_i, delim_whitespace=True, header=0, index_col=0)

                else:
                    print('-- Warning NAN flux entries')

                counter += 1

# Loop through the dicts and get the data:
lineMeanSampleProperties('O3_5007A', df_dict)
lineMeanSampleProperties('O3_4959A', df_dict)
lineMeanSampleProperties('S3_6312A', df_dict)
lineMeanSampleProperties('He1_5016A', df_dict)
lineMeanSampleProperties('H1_8392A', df_dict)

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