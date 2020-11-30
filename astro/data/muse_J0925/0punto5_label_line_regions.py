import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams
from astropy.wcs import WCS
from src.specsiser.print.plot import STANDARD_PLOT
from astro.data.muse.common_methods import compute_line_flux_image, image_array_binning

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

for i, obj in enumerate(objList):

    # Data location
    cube_address_i = dataFolder/fileList[i]
    objFolder = resultsFolder
    db_addresss = resultsFolder/f'{obj}_database.txt'

    # Load data
    wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')
    wave = wave / (1 + z_objs[i])
    obj_db = pd.read_csv(db_addresss, delim_whitespace=True, header=0, index_col=0)
    print(f'\n- {obj}: Cube dimensions {cube.shape}')

    # Get line region data
    lineFlux_dict, levelFlux_dict, levelText_dict = compute_line_flux_image(lineAreas,
                                                                            cube,
                                                                            z_objs[i],
                                                                            percent_array=pertil_array)

    # Reconstruct the max cube
    max_cube = cube[:, :, :].data.data
    image_max = np.nanmax(max_cube, axis=0)

    # Reconstruct the max flux from file
    cube_shape = obsData['sample_data']['cube_size_array']
    txt_max = np.reshape(obj_db['max_flux'].values, cube_shape.astype(int))
    print(image_max.shape, txt_max.shape)

    idx_pair = (45, 40)
    idx_j, idx_i = idx_pair
    idx_db = (obj_db.y_voxel == idx_j) & (obj_db.x_voxel == idx_i)
    # print(np.unravel_index(counter, matrix.shape))

    print(obj_db.loc[idx_db])
    print(image_max[idx_j, idx_i], txt_max[idx_j, idx_i], obj_db.loc[idx_db, 'max_flux'].values)
    idx_true = obj_db.max_flux == 99.190002
    index_num = obj_db.loc[idx_true].index.values[0]
    print(obj_db.loc[idx_true, 'max_flux'])