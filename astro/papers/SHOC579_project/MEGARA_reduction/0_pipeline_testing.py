import shutil

import lime
import numpy as np
from astropy.io import fits
from pathlib import Path
import lime as lm
import os
from shutil import copyfile
import pandas as pd

# Loading project configuration file
obs_conf = lm.load_cfg(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\obsConf.ini')
reduction_cfg = obs_conf['Megara_reduction']

# Dataframe with files list
rd_df_address = Path(reduction_cfg['rd_df_address'])
# sample_DF = lm.load_lines_log(f'{rd_df_address}.txt')

# Stating folder structure
instructions_folder = Path(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\MEGARA_reduction')
reduction_folder = Path(reduction_cfg['root_folder'])
data_folder = reduction_folder/'data'
source_folder = f'S:\Astro_data\Observations\SHOC579\MEGARA\orig_files\GTC27-21B'

# Checking the files
sample_DF = pd.DataFrame(columns=['OB', 'index_type', 'object', 'type', 'reduc_tag', 'VPH', 'address', 'old_address'])
rename_dict = {}

for root, dirs, files in os.walk(source_folder):
    for i, old_name in enumerate(files):
        if old_name.endswith('.fits'):

            # Renaming data
            OB_name = root[root.find('OB'):root.find('OB')+6]
            type_name = root[root.rfind("\\")+1:]

            with fits.open(f'{root}\{old_name}') as hdu_list:
                data = hdu_list[0].data
                hdr = hdu_list[0].header

            VPH = fits.getval(f'{root}\{old_name}', 'VPH')
            obj = fits.getval(f'{root}\{old_name}', 'OBJECT')

            # Moving the file
            if obj == 'SHOC579':
                ref_obj = 'SHOC579'
            elif obj == 'SPSTD_HR8634':
                ref_obj = 'HR8634'
            elif obj == 'HR7596':
                ref_obj = 'HR7596'
            else:
                ref_obj = type_name

            new_name = f'{OB_name}_{ref_obj}_{i}.fits'
            rename_dict[f'{root}\{old_name}'] = new_name

            print(f'{root}\{old_name}', f'->', f'{data_folder}\{new_name}')
            shutil.copyfile(f'{root}\{old_name}', f'{data_folder}\{new_name}')

            # Saving the files data to a dataframe
            sample_DF.loc[new_name, 'OB'] = OB_name
            sample_DF.loc[new_name, 'index_type'] = i
            sample_DF.loc[new_name, 'object'] = ref_obj
            sample_DF.loc[new_name, 'VPH'] = VPH
            sample_DF.loc[new_name, 'type'] = type_name
            sample_DF.loc[new_name, 'reduc_tag'] = 'raw'
            sample_DF.loc[new_name, 'address'] = f'{data_folder}\{new_name}'
            sample_DF.loc[new_name, 'old_address'] = f'{root}\{old_name}'

lime.save_line_log(sample_DF, rd_df_address)