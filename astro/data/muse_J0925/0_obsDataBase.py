import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
import astropy.units as u
from mpdaf.obj import deg2sexa

lineAreas = {'H1_6563A': (6533.0, 6596.0),
             'S3_6312A': (6310.0, 6319.0),
             'O3_5007A': (4999.0, 5025.0),
             'O3_4363A': (4355.0, 4374.0)}

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
# db_format = {'DEC_deg': lambda x: '%.8f' % x, 'RA_deg': lambda x: '%.8f' % x}

for idx, obj in enumerate(objList):

    # Specify inputs and output
    cube_address = dataFolder/fileList[idx]
    objFolder = resultsFolder
    db_addresss = objFolder/f'{obj}_database.txt'

    # Load the data
    wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')
    wave = wave / (1 + z_objs[idx])
    print(f'\n- {obj}: Cube dimensions {cube.shape}')

    # Object database as a pandas dataframe
    size_y, size_x = cube.shape[1], cube.shape[2]
    obsDF = pd.DataFrame(index=np.arange(size_y*size_x), columns=db_headers)

    # Voxel indeces where Cube shape = (lambda , Y, X)
    y_range, x_range = np.ogrid[0:size_y, 0:size_x]
    X, Y = np.meshgrid(x_range, y_range)
    X_flatten, Y_flatten = X.flatten(), Y.flatten()

    voxel_label_array = np.empty(obsDF.index.size, dtype=str)
    voxel_label_array.fill('-')
    voxel_label_array = np.core.defchararray.add(X_flatten.astype(str), voxel_label_array)
    voxel_label_array = np.core.defchararray.add(voxel_label_array, Y_flatten.astype(str))

    obsDF['y_voxel'] = Y_flatten
    obsDF['x_voxel'] = X_flatten
    obsDF['voxel_label'] = voxel_label_array

    obsDF['max_flux'] = np.nan
    max_cube = cube[:, :, :].data.data
    max_matrix = np.nanmax(max_cube, axis=0)
    obsDF['max_flux'] = max_matrix.flatten()

    # voxel_label_array = np.empty(obsDF.index.size, dtype=str)
    # voxel_label_array.fill('-')
    # voxel_label_array = np.core.defchararray.add(Y_flatten.astype(str), voxel_label_array)
    # voxel_label_array = np.core.defchararray.add(voxel_label_array, X_flatten.astype(str))
    # obsDF['voxel_label'] = voxel_label_array

    # Voxel coordinates
    DEC_0, RA_0 = cube.wcs.get_start()
    DEC_1, RA_1 = cube.wcs.get_end()
    DEC_arange = np.linspace(DEC_0, DEC_1, cube.shape[2])
    RA_arange = np.linspace(RA_0, RA_1, cube.shape[1])
    coord_array = np.c_[Y_flatten, X_flatten]
    deg_array = cube.wcs.pix2sky(coord_array, unit=u.deg)
    sexa_array = deg2sexa(deg_array) # This operation takes some time
    obsDF['DEC'] = sexa_array[:, 0]
    obsDF['RA'] = sexa_array[:, 1]
    obsDF['DEC_deg'] = deg_array[:, 0]
    obsDF['RA_deg'] = deg_array[:, 1]

    # Save the dataframe
    with open(db_addresss, 'wb') as output_file:
        string_DF = obsDF.to_string()
        output_file.write(string_DF.encode('UTF-8'))
