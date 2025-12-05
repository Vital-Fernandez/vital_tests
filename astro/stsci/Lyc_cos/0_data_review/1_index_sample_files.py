from pathlib import Path
import pandas as pd
import lime
from astro.stsci.tools import list_files_with_extension, add_cos_obs


# Data location
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
obs_folder = Path('/home/vital/Astrodata/STScI')
LyC_obs_data = obs_folder/'LyC_leakers_COS'/'Direct_downloads'
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

# Empty dataframe
hdrs = ['sample', 'id', 'offset_id', 'state', 'object', 'PID', 'RA', 'DEC', 'redshift',
        'instr', 'grating', 'cenwave', 'life_adj', 'distap', 'detector',
        'filecode', 'filepath']

sample_df = pd.DataFrame(columns=hdrs)

# Calibrated files
ext = '_x1dsum.fits'
list_x1dsum = list_files_with_extension(LyC_obs_data, ext)
add_cos_obs(list_x1dsum, sample_df, 'LyC_cos', 'x1dsum', cfg_sample)

# HASP files
ext = '_cspec.fits'
list_hasp = list_files_with_extension(LyC_obs_data, ext)
add_cos_obs(list_hasp, sample_df, 'LyC_cos', 'hasp_mast', cfg_sample)

# Manual HASP files
ext = '_x1d.fits'
list_files = list_files_with_extension(LyC_obs_data, ext)
add_cos_obs(list_files, sample_df, 'LyC_cos', 'x1d', cfg_sample)

# Manual rawacq files
ext = '_flt.fits'
list_files = list_files_with_extension(LyC_obs_data, ext)
add_cos_obs(list_files, sample_df, 'LyC_cos', 'flt', cfg_sample)

# Manual cal files
ext = '_cal.fits'
list_files = list_files_with_extension(LyC_obs_data, ext)
add_cos_obs(list_files, sample_df, 'LyC_cos', 'cal', cfg_sample)

# Manual cal files
ext = '_mos.fits'
list_files = list_files_with_extension(LyC_obs_data, ext)
add_cos_obs(list_files, sample_df, 'LyC_cos', 'mos', cfg_sample)

# Manual cal files
ext = '_drz.fits'
list_files = list_files_with_extension(LyC_obs_data, ext)
add_cos_obs(list_files, sample_df, 'LyC_cos', 'drz', cfg_sample)

print(f'Duplicated identifiers')
groups = sample_df.groupby(list(sample_df.columns[:4])).groups
duplicate_groups = {k: list(v) for k, v in groups.items() if len(v) > 1}
print(duplicate_groups)
print()

# Save the dataframe
lime.save_frame(project_folder/'stsci_samples_v0.csv', sample_df)
lime.save_frame(project_folder/'stsci_samples_v0.txt', sample_df)