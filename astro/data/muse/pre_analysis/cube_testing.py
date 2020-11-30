import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
import astropy.units as u
from mpdaf.obj import deg2sexa

# Declare data and files location
conf_file_address = '../muse_greenpeas.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)
objList = obsData['sample_data']['object_list']
fileList = obsData['sample_data']['file_list']
dataFolder = Path(obsData['sample_data']['data_folder'])
results_folder = Path(obsData['sample_data']['results_folder'])
z_objs = obsData['sample_data']['z_array']
db_headers = obsData['sample_data']['database_header_list']
db_format = {'DEC_deg': '{: 0.8f}', 'RA_deg': '{: 0.8f}'}
# db_format = {'DEC_deg': lambda x: '%.8f' % x, 'RA_deg': lambda x: '%.8f' % x}

for idx, obj in enumerate(objList):

    # Specify inputs and output
    cube_address = dataFolder/fileList[idx]
    objFolder = results_folder/f'{obj}'
    db_addresss = objFolder/f'{obj}_database.txt'

    # Load the data
    wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')
    wave = wave / (1 + z_objs[idx])
    print(f'\n- {obj}: Cube dimensions {cube.shape}')