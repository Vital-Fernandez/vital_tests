from pathlib import Path
import pandas as pd
import lime
from lime.fitting.lines import c_KMpS
import re
import numpy as np
import os

'''
This script reads the FWHM vs wavelength values estimated from
the COS LSF, and calculates the total observed FWHM:

FWHM_spec^2 = FWHM_LSF^2 + (FWHM_img^2 - FWHM_PSF^2)

where 
FWHM_img = FWHM of the image profile measured in the NUV TA frames
FWHM_PSF = FWHM of the PSF in the NUV TA images ~2 pixels
'''

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
lsf_folder =  Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_file')
lsf_fittings_path = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_fitted')
fwhm_total_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_total')
project_folder = Path('/')

# Cfg file
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v0.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filecode"].str.contains(pattern, na=False)]

# Get the targe
obj_fwhm_dict = cfg_sample['LyC_acq_image_fwhm_pixels']
targ_arr = np.array(obj_fwhm_dict.keys())
fwhm_img_arr =  np.array(obj_fwhm_dict.values())

# Analysis configuraiton
FWHM_PSF = 2 # NUV detector
disp_fator = {"G130M": 9.93e-3, "G160M": 12.23e-3, "nuv": 37.0e-3}
grat_arr = np.array(['G130M', 'G160M', 'G185M'])

# Loop throught the objects
for i, (targ, fwhm_img) in enumerate(obj_fwhm_dict.items()):

    # Loop through the gratings
    fls_file_list = []
    for grat in grat_arr:
        if grat != 'G185M':
            if ('Haro11' in targ) or ('Izw18' in targ):
                sub_labels = cfg_sample['multi_target_labels'][targ.split('_')[0]][targ]
                idcs = sample_df.index.get_level_values('id').isin(sub_labels)
                idcs = idcs & (sample_df.grating == grat) & pd.notnull(sample_df.life_adj)
            else:
                idcs = (sample_df.object == targ) & (sample_df.grating == grat) & pd.notnull(sample_df.life_adj)
            obj_df = sample_df.loc[idcs, ['grating', 'cenwave', 'life_adj']]
            uniq_comb = obj_df.drop_duplicates()
            for idx in uniq_comb.index:
                lsf_name = f'aa_LSFTable_{grat}_{uniq_comb.loc[idx].cenwave}_LP{int(uniq_comb.loc[idx].life_adj)}_cn.dat'
                fls_file_list.append(lsf_name)
        else:
            if 'nuv_model_lsf.dat' not in fls_file_list:
                fls_file_list.append('nuv_model_lsf.dat')

    # Loop through the matched LSF
    for lsf_path in fls_file_list:
        wave, fwhm_lsf = np.loadtxt(lsf_fittings_path/f'{Path(lsf_path).stem}_fitted.dat', usecols=(0, 1), unpack=True)
        grating =  lsf_path.split('_')[2] if 'nuv' not in lsf_path else 'nuv'
        cenwave = lsf_path.split('_')[3] if 'nuv' not in lsf_path else None
        life_apj = lsf_path.split('_')[4] if 'nuv' not in lsf_path else None

        # --- allocate tables (first row = wavelength) ---
        table_pixels = np.zeros((2, wave.size))
        table_pixels[0, :] = wave

        table_vel = np.zeros((2, wave.size))
        table_vel[0, :] = wave

        # Total FWHM
        fwhm_spec = np.sqrt(fwhm_lsf**2 + (fwhm_img**2 - FWHM_PSF**2))

        # pixels -> Angstroms
        fwhm_spec_A = fwhm_spec * disp_fator.get(grating, 1)

        # Angstroms -> velocity
        fwhm_spec_vel = (fwhm_spec_A * c_KMpS) / wave

        table_pixels[1, :] = fwhm_spec
        table_vel[1, :] = fwhm_spec_vel

        # --- save to ASCII (dynamic columns & names) ---
        # cols_pixels = [table_pixels[0, :]] + [table_pixels[i + 1, :] for i in range(len(targ))]
        # cols_vel = [table_vel[0, :]] + [table_vel[1, :] for i in range(len(targ))]
        # names = ['#wave'] + list(targ)

        # Build filenames near the input file; tweak to your preferred pattern
        base = f'{targ}_{grating}_{cenwave}_{life_apj}'  # e.g. G130M_1291_LP3_fwhm
        out_pixels = os.path.join(os.path.dirname(lsf_path), f"{base}_total_pixels.dat")
        out_vel = os.path.join(os.path.dirname(lsf_path), f"{base}_total_velocity_kms.dat")

        np.savetxt(fwhm_total_folder/f"{base}_pixel.txt", table_pixels.T, fmt="%.6f", delimiter=" ", header=f"wave {targ}", comments="")
        np.savetxt(fwhm_total_folder/f"{base}_vel.txt", table_vel.T, fmt="%.6f", delimiter=" ", header=f"wave {targ}", comments="")

        print(f"Wrote: {out_pixels}")
        print(f"Wrote: {out_vel}")