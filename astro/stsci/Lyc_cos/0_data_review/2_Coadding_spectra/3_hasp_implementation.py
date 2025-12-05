import re
import lime
import shutil
from pathlib import Path
from astro.stsci.tools import list_files_with_extension, add_cos_obs, run_hasp_wrapper, move_files


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v0.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filecode"].str.contains(pattern, na=False)]
sample_df = sample_df.loc[sample_df.index.get_level_values('state') == 'x1d']

# Get list of targets
target_list = sample_df.object.unique()
grating_list = ['G130M', 'G160M', 'G185M']

# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    # Object folders the files
    input_folder_single = obs_folder/'LyC_leakers_COS'/'objects_x1d'/f'{obj}'
    output_folder_single = obs_folder/'LyC_leakers_COS'/'obj_hasp'/f'{obj}'

    # # Move the files
    # idcs_x1d = (sample_df.index.get_level_values('state') == 'x1d') & (sample_df.object == obj) & (sample_df.grating.isin(grating_list))
    # file_arr = sample_df.loc[idcs_x1d, 'filepath'].to_numpy()
    # move_files(file_arr, obs_folder, input_folder_single)
    #
    # # Run HASP wrapper
    # if output_folder_single.exists():
    #     shutil.rmtree(output_folder_single)
    # output_folder_single.mkdir(parents=True, exist_ok=True)
    # run_hasp_wrapper(input_folder_single, output_folder_single, cross_program=True)

    # Append new files to table
    ext = '_cspec.fits'
    list_files = list_files_with_extension(output_folder_single, ext)
    add_cos_obs(list_files, sample_df, 'LyC_cos', 'cspec_manual', cfg_sample)

    ext = '_aspec.fits'
    list_files = list_files_with_extension(output_folder_single, ext)
    add_cos_obs(list_files, sample_df, 'LyC_cos', 'aspec_manual', cfg_sample)

# Save the dataframe
lime.save_frame(project_folder / 'stsci_samples_v1.csv', sample_df)
lime.save_frame(project_folder / 'stsci_samples_v1.txt', sample_df)





# # Multi PID objects
# obj_list_multID = group_PID[group_PID > 1].index.tolist()
# for i, obj_name in enumerate(obj_list_multID):
#     idcs_x1d = (sample_df.index.get_level_values('state') == 'x1d') & (sample_df.object == obj_name)
#     file_arr = sample_df.loc[idcs_x1d, 'filepath'].to_numpy()
#     move_files(file_arr, obs_folder, obj_folder_list[i])

# # Run wrapper for multiple PI
# if output_folder_mult.exists():
#     shutil.rmtree(output_folder_mult)
# output_folder_mult.mkdir(parents=True, exist_ok=True)
# for obj_folder in obj_folder_list:
#     run_hasp_wrapper(obj_folder, output_folder_mult, cross_program=True)

# list_hasp = list_files_with_extension(output_folder_mult, '_cspec.fits')
# add_cos_obs(list_hasp, sample_df, 'LyC_cos', 'hasp_mult', cfg_sample)
# lime.save_frame(project_folder/'stsci_samples_v1.csv', sample_df)

# list_hasp = list_files_with_extension(output_folder_mult, '_aspec.fits')
# add_cos_obs(list_hasp, sample_df, 'LyC_cos', 'hasp_aspec', cfg_sample)
# lime.save_frame(project_folder/'stsci_samples_v1.csv', sample_df)