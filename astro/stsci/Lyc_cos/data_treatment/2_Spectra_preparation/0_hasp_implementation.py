import re

import matplotlib.pyplot as plt

import lime
import shutil
from pathlib import Path
from astro.stsci.tools import list_files_with_extension, add_cos_obs, run_hasp_wrapper, move_files
import lineid_plot


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v0.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Get list of targets
target_list = sample_df.object.unique()
grating_list = ['G130M', 'G160M', 'G185M']

# Lines list
line_labels = cfg_sample['line_properties']['obs_lines']
particles, trans_arr, latex_arr = lime.label_decomposition(line_labels)

# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    print(f'{i}) Galaxy {obj}')
    input_folder_single = obs_folder/'LyC_leakers_COS'/'objects_x1d'/f'{obj}'
    output_folder_single = obs_folder/'LyC_leakers_COS'/'obj_hasp'/f'{obj}'

    # Move the files
    idcs_x1d = (sample_df.index.get_level_values('state') == 'x1d') & (sample_df.object == obj) & (sample_df.grating.isin(grating_list))
    file_arr = sample_df.loc[idcs_x1d, 'filepath'].to_numpy()
    # move_files(file_arr, obs_folder, input_folder_single)

    # # Run HASP wrapper
    # if output_folder_single.exists():
    #     shutil.rmtree(output_folder_single)
    # output_folder_single.mkdir(parents=True, exist_ok=True)
    # run_hasp_wrapper(input_folder_single, output_folder_single, cross_program=True)

    # Append new files to table
    ext = '_cspec.fits'
    list_files = list_files_with_extension(output_folder_single, ext)
    add_cos_obs(list_files, sample_df, 'LyC_cos', 'hasp_manual', cfg_sample)

    ext = '_aspec.fits'
    list_files = list_files_with_extension(output_folder_single, ext)
    add_cos_obs(list_files, sample_df, 'LyC_cos', 'aspec_manual', cfg_sample)

# Save the dataframe
lime.save_frame(project_folder / 'stsci_samples_v1.csv', sample_df)
lime.save_frame(project_folder / 'stsci_samples_v1.txt', sample_df)