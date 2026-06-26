from pathlib import Path
import pandas as pd
import lime
from astro.stsci.tools import list_files_with_extension, add_cos_obs


# Data location
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
obs_folder = Path('/home/vital/Astrodata/STScI')
LyC_obs_data = obs_folder/'LyC_leakers_COS'/'Direct_downloads'
cfg_sample = lime.load_cfg('../Lyc_cos.toml')

# Empty dataframe
hdrs = ['sample', 'id', 'offset_id', 'state', 'object', 'PID', 'RA', 'DEC', 'redshift',
        'instr', 'grating', 'cenwave', 'life_adj', 'distap', 'detector',
        'filecode', 'filepath']

sample_df = pd.DataFrame(columns=hdrs)

# Calibrated files
ext = '_aspec.fits'
fits_folder = '/home/vital/Astrodata/STScI/LyC_leakers_COS/obj_hasp'
list_files = list_files_with_extension(fits_folder, ext)
add_cos_obs(list_files, sample_df, 'LyC_cos', 'aspec', cfg_sample)

print(f'Duplicated identifiers')
groups = sample_df.groupby(list(sample_df.columns[:4])).groups
duplicate_groups = {k: list(v) for k, v in groups.items() if len(v) > 1}
print(duplicate_groups)
print()

# Add the columns with additional spectra files
sample_df['LyAlpha_fitting'] = 'none'
sample_df['metals_fitting'] = 'none'
for index in sample_df['object'].index:
        sample_df.loc[index, 'filepath'] = f'LyC_leakers_COS/{Path(sample_df.loc[index, 'filepath']).name}'
        sample_df.loc[index, 'LyAlpha_fitting'] = f'LyC_leakers_COS/{sample_df.loc[index].object}_LyAlpha_spec.txt'
        sample_df.loc[index, 'metals_fitting'] = f'LyC_leakers_COS/{sample_df.loc[index].object}_metals_spec.txt'

# Save the dataframe
lime.save_frame(project_folder/'stsci_LyC_HST_specsy_v0.csv', sample_df)
lime.save_frame(project_folder/'stsci_LyC_HST_specsy_v0.txt', sample_df)
print(f'Saving at {project_folder/'stsci_LyC_HST_specsy_v0.csv'}')
print('Unique objects')
target_lits = sample_df.object.unique()
print(target_lits)
print(f'{len(target_lits)} targets')
print(f'Repeated entries: {sample_df[sample_df["object"].duplicated(keep=False)].index.to_numpy()}')