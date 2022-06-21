import os
import shutil
import pandas as pd
import lime as lm
from astropy.io import fits
from pathlib import Path

# Configuration file
cfg_file = '../obsConf.ini'
obs_conf = lm.load_cfg(Path(cfg_file))
reduction_cfg = obs_conf['Megara_reduction']

# Data location
reduction_folder = Path(reduction_cfg['root_folder'])
data_folder = reduction_folder/'data'
rd_df_address = Path(reduction_cfg['rd_df_address'])
source_folder = Path(reduction_cfg['fits_folder'])

# Checking the files
sample_DF = pd.DataFrame(columns=['OB', 'index_type', 'object', 'type', 'reduc_tag', 'VPH', 'address', 'old_address'])
rename_dict = {}

files_counter = 0
for root, dirs, files in os.walk(source_folder):
    for i, old_name in enumerate(files):
        if old_name.endswith('.fits'):

            # Renaming data
            OB_name = root[root.find('OB'):root.find('OB')+6]
            type_name = root[root.rfind("/")+1:]

            with fits.open(f'{root}/{old_name}') as hdu_list:
                data = hdu_list[0].data
                hdr = hdu_list[0].header

            VPH = fits.getval(f'{root}/{old_name}', 'VPH')
            obj = fits.getval(f'{root}/{old_name}', 'OBJECT')

            # Moving the file
            if obj == 'SHOC579':            # Target
                ref_obj = 'SHOC579'
            elif obj == 'SPSTD_HR8634':     # Standard star
                ref_obj = 'HR8634'
            elif obj == 'HR7596':           # Standard star
                ref_obj = 'HR7596'
            else:
                ref_obj = type_name

            new_name = f'{OB_name}_{VPH}_{ref_obj}_{i}.fits'
            rename_dict[f'{root}/{old_name}'] = new_name

            # Copy the files to another folder with a new name
            old_path = Path(root)/old_name
            new_path = data_folder/new_name
            files_counter += 1
            print(f'{files_counter}) {old_name}', f'->', f'{new_path} {"" if new_path.is_file() else "(False)"}')
            # shutil.copyfile(f'{old_path}', f'{new_path}')

            # Saving the files data to a dataframe
            sample_DF.loc[new_name, 'OB'] = OB_name if OB_name != 'OB0004' else 'OB0003'
            sample_DF.loc[new_name, 'index_type'] = i
            sample_DF.loc[new_name, 'object'] = ref_obj
            sample_DF.loc[new_name, 'VPH'] = VPH
            sample_DF.loc[new_name, 'type'] = type_name
            sample_DF.loc[new_name, 'reduc_tag'] = 'raw'
            sample_DF.loc[new_name, 'address'] = f'{new_path}'
            sample_DF.loc[new_name, 'old_address'] = f'{old_path}'

print(rd_df_address)
sample_DF.sort_values(['OB'], ascending=True, inplace=True)
lm.save_line_log(sample_DF, rd_df_address)
