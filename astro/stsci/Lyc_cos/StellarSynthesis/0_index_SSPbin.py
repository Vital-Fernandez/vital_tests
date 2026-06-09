from pathlib import Path

import lime
from pandas import DataFrame

hdrs = ['z', 'age', 'imf', 'Mmin', 'Mmax',  'fname']
bin_df = DataFrame(columns=hdrs)

ssp_folder = Path('/home/vital/Astrodata/BPASS_v2.3/spectra-bin_byrne23')
out_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs')

# Loop through the files
bin_list = list(ssp_folder.glob('*.dat'))
for fpath in bin_list:
    fname = fpath.stem
    parts = fname.split('.')

    m_low, m_high = parts[0].split('-')[2].split('_')
    m_low = m_low[3:]

    imf = parts[1]

    metal = f'{parts[2][3:]}'if parts[2].startswith('zem') else float(f'0.{parts[2][1:]}')

    age = parts[3][3:]

    print(fname, '  ->  ', imf, metal, age)
    bin_df.loc[len(bin_df)] = {'z': metal, 'age': age, 'imf': imf, 'Mmin': m_low, 'Mmax': m_high, 'fname': fname}

# Sort
bin_df = bin_df.sort_values(['imf', 'z', 'age']).reset_index(drop=True)

# Save
dfname =out_folder/'BPASS_v2.3.1.txt'
lime.save_frame(dfname, bin_df)

