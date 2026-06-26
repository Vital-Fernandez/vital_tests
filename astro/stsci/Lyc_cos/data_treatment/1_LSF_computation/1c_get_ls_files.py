from pathlib import Path
import pandas as pd
import lime
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np
import pdfkit
from astro.stsci.tools import fetch_files
import re

# Data location
project_folder = Path('/')
obs_folder = Path('/home/vital/Astrodata/STScI')
lsf_folder =  Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_file')
# https://spacetelescope.github.io/hst_notebooks/notebooks/COS/LSF/LSF.html#undL

# Cfg file
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v0.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filecode"].str.contains(pattern, na=False)]

# Get list of x1d files:
idcs_x1dm = sample_df.index.get_level_values('state') == 'x1d'
df_hdrs = pd.DataFrame(columns=["DETECTOR", "OPT_ELEM", "LIFE_ADJ", "CENWAVE"])
for i, idx in enumerate(sample_df.loc[idcs_x1dm].index):

    file_path = obs_folder/sample_df.loc[idx, 'filepath']
    hdr = fits.getheader(file_path, ext=0)

    keywords = ["DETECTOR", "OPT_ELEM", "LIFE_ADJ", "CENWAVE", "DISPTAB"]
    df_hdrs.loc[i,'DETECTOR'] = hdr['DETECTOR']
    df_hdrs.loc[i,'OPT_ELEM'] = hdr['OPT_ELEM']
    df_hdrs.loc[i,'LIFE_ADJ'] = hdr['LIFE_ADJ']
    df_hdrs.loc[i,'CENWAVE'] = hdr['CENWAVE']

    fetch_files(det=hdr['DETECTOR'], grating=hdr['OPT_ELEM'], lpPos=hdr['LIFE_ADJ'], cenwave=hdr['CENWAVE'],
                disptab=hdr['DISPTAB'], datadir=lsf_folder)

# print('')
# unique_rows = df_hdrs[~df_hdrs.duplicated(keep=False)]
# n_unique_rows = len(unique_rows)
# print(f'Unique rows: {n_unique_rows} of {df_hdrs.index.size}')
#

