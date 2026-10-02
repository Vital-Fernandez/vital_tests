import re
import lime
import numpy as np
from pathlib import Path

lime.theme.set_style('dark')


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Order for objects
target_list = sample_df.object.unique()
target_list = ['SBS0335052', 'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
               'UGCA281', 'NGC1705',  'NGC4861',  'Haro2', 'He2-10', 'MRK1450',         # Codo
               'UGC4483', 'VIIZw403',  'NGC2366',                                       # Dent
               'Haro11_A', 'Haro11_B', 'Haro11_C']

stellar_lines = cfg_sample['default_line_fitting']['stellar_lines']

# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    if i == 3:

        # Object folders the files
        spec_folder = project_folder/'LyC_leakers_COS'/'metals_line_fitting'
        stellar_folder = project_folder/'LyC_leakers_COS'/'SSPs'
        print(f'{i}) Galaxy {obj}')

        fname = spec_folder/f'{obj}_rebinned_LymanAlphaCorrected.txt'
        spec = lime.Spectrum.from_file(fname, instrument='text')

        fname = stellar_folder/'stellar_masks'/f'{obj}_mask_intervals.txt'
        spec.check.masks(fname, intvls= stellar_folder/'stellar_masks'/'template_mask_interval.txt', line_list=stellar_lines,
                         ax_cfg={'title':obj}, rest_frame=True)
        # spec.plot.spectrum(line_list=stellar_lines)
