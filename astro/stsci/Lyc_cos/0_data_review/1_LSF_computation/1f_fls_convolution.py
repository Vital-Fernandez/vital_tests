from pathlib import Path
import numpy as np
import pandas as pd
import lime
import re

# Data location
project_folder = Path('/')
obs_folder = Path('/home/vital/Astrodata/STScI')
lsf_folder =  Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_file')
lsf_fittings_path = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_fitted')
fwhm_total_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_total')
convolved_lsf_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_convolved')

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

# Loop throught the objects
for i, (targ, fwhm_img) in enumerate(obj_fwhm_dict.items()):

    if ('Haro11' in targ) or ('Izw18' in targ):
        sub_labels = cfg_sample['multi_target_labels'][targ.split('_')[0]][targ]
        idcs = sample_df.index.get_level_values('id').isin(sub_labels)
    else:
        idcs = (sample_df.object == targ)

    # Get object observations
    idcs = idcs & (sample_df.index.get_level_values('state') == 'x1d') & pd.notnull(sample_df.life_adj)
    obj_df = sample_df.loc[idcs, ['grating', 'cenwave', 'life_adj']]

    # Loop throught the observation combinations
    uniq_comb = obj_df.drop_duplicates()
    for idx in uniq_comb.index:
        grating, cenwave, life_adj = uniq_comb.loc[idx].values

        # Open the COS LSF
        lsf_file = f'aa_LSFTable_{grating}_{cenwave}_LP{int(life_adj)}_cn.dat' if grating != 'G185M' else 'nuv_model_lsf.dat'
        matrix_arr = np.loadtxt(lsf_folder/lsf_file)
        wave_lsf, lsf_matrix = matrix_arr[0, :], matrix_arr[1:, :]

        # Open the fitted LSF function
        lsf_fitted_file = f'aa_LSFTable_{grating}_{cenwave}_LP{int(life_adj)}_cn_fitted.dat' if grating != 'G185M' else 'nuv_model_lsf_fitted.dat'
        wave_cos, fwhm_cos = np.loadtxt(lsf_fittings_path/lsf_fitted_file, unpack=True)

        # Open the total object total FWHM
        fwhmTotal_pixel_fname = f'{targ}_{grating}_{cenwave}_LP{int(life_adj)}_pixel.txt' if grating != 'G185M' else f'{targ}_nuv_None_None_pixel.txt'
        wave_fwhm, fwhm_object = np.loadtxt(fwhm_total_folder/fwhmTotal_pixel_fname, unpack=True, skiprows=1)

        # Interpolate to the original lsf wavelength
        object_fwhm_interp = np.interp(wave_lsf, wave_fwhm, fwhm_object)
        cos_fwhm_interp = np.interp(wave_lsf, wave_cos, fwhm_cos)

        # Loop throught the wavelength values
        new_lsf = np.zeros(lsf_matrix.shape)
        sigma_corr_arr = np.sqrt((object_fwhm_interp / 2.35482) ** 2 - (cos_fwhm_interp / 2.35482) ** 2)
        # for j, wavelength in enumerate(wave_lsf):
        #
        #     ### Convolve COS LSF with FWHM_total ####
        #     kernel = Gaussian1DKernel(stddev=sigma_corr_arr[j])
        #     new_lsf[:, j] = convolve(lsf_matrix[:, j], kernel)
        #
        #     fig, ax = plt.subplots()
        #     ax.plot(lsf_matrix[:, j], label=f'Original')
        #     ax.plot(new_lsf[:, j], label=f'Convolved')
        #     ax.update({'title': f'{targ} LSF convolution for wavelength {wavelength}Å \n ACQ image FWHM = {fwhm_img}, {grating}, cenwave = {cenwave},'
        #                         f' LP{int(life_adj)}'})
        #     ax.legend()
        #     plt.show()

        # Save to a text file
        lsf_convolved_file = f'{targ}_aa_LSFTable_{grating}_{cenwave}_LP{int(life_adj)}_cn_convolved.dat' if grating != 'G185M' else f'{targ}_nuv_model_lsf_convolved.dat'
        np.savetxt(convolved_lsf_folder/lsf_convolved_file, new_lsf)