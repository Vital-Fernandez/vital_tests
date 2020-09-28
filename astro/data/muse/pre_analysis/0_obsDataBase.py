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
db_headers = obsData['sample_data']['database_header_list']
# db_format = {'DEC_deg': '{: 2.5f}', 'RA_deg': '{: 2.5f}'}
db_format = {'DEC_deg': lambda x: '%.8f' % x, 'RA_deg': lambda x: '%.8f' % x}


# formats = {'c': '{: 10d}', 'd': '{: 2.5f}'}


for i, obj in enumerate(objList):

    # Load the data
    print(f'\n- {obj}')
    file_address_i = f'{dataFolder}/{fileList[i]}'
    wave, cube, header = sr.import_fits_data(file_address_i, instrument='MUSE')
    wave = wave / (1 + z_objs[i])

    # Cube shape(lambda , Y, X)
    size_y, size_x = cube.shape[1], cube.shape[2]
    obsDF = pd.DataFrame(index=np.arange(size_y*size_x), columns=db_headers)

    idx_counter = 0
    for idx_j in np.arange(size_y):
        for idx_i in np.arange(size_x):

            # Voxel references
            voxel_label = f'{idx_j}_{idx_i}'
            obsDF.loc[idx_counter, 'voxel_label'] = voxel_label
            obsDF.loc[idx_counter, 'y_voxel'] = idx_j
            obsDF.loc[idx_counter, 'x_voxel'] = idx_i

            # Voxel coordinates
            coord_sky = cube.wcs.pix2sky((idx_j, idx_i), unit=u.deg)
            dec, ra = deg2sexa(coord_sky)[0]
            obsDF.loc[idx_counter, 'DEC'] = dec
            obsDF.loc[idx_counter, 'RA'] = ra
            obsDF.loc[idx_counter, 'DEC_deg'] = coord_sky[0][0]
            obsDF.loc[idx_counter, 'RA_deg'] = coord_sky[0][1]
            print(f'-- Voxel {voxel_label}: {coord_sky[0][0]} {coord_sky[0][1]}')

            # Increase voxel counter reference
            idx_counter += 1

    # Save the dataframe
    table_name = fileList[i].replace('.fits', '_database.txt')
    table_address = f'{dataFolder}/{table_name}'
    with open(table_address, 'wb') as output_file:
        string_DF = obsDF.to_string(formatters=db_format)
        output_file.write(string_DF.encode('UTF-8'))

    # # Get astronomical coordinates one pixel
    # coord_sky = cube.wcs.pix2sky(idx_voxel, unit=u.deg)
    # dec, ra = deg2sexa(coord_sky)[0]
    # wcs_cube = WCS(cube.data_header)