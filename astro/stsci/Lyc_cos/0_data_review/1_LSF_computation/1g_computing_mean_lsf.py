import re
import numpy as np
import lime
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib.widgets import SpanSelector
from lime.fitting.lines import c_KMpS
from astro.stsci.tools import on_click, IntervalSelector
from astropy.io import fits
from scipy.interpolate import interp1d

lime.theme.set_style('dark')

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
output_folder = project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'
lsf_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/lsf_convolved')
lsf_obj_folder = project_folder/'LyC_leakers_COS'/'lsf_objects'

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Order for objects
target_list = sample_df.object.unique()
target_list = cfg_sample['galaxy_sample']['sn_sort_list']

# Loop through the targets and generate the HASP files
check_masks = False
for i, obj in enumerate(target_list):

    if i >= 15:

        # Object folders the files
        print(f'{i}) Galaxy {obj}, z = {cfg_sample['Galaxy_redshifts'][obj]}')
        input_folder_single = obs_folder / 'LyC_leakers_COS' / 'objects_x1d' / f'{obj}'
        output_folder_single = obs_folder / 'LyC_leakers_COS' / 'obj_hasp' / f'{obj}'

        # Index the files
        idcs_out = (sample_df.object == obj) & (sample_df.index.get_level_values('state') == 'aspec_manual')

        # Load the files
        pname = obs_folder / sample_df.loc[idcs_out].filepath[0]
        spec = lime.Spectrum.from_file(pname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj], norm_flux=1e-17)
        data = fits.getdata(pname, ext=2)

        # Rebinned
        spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)
        wave_interp = np.arange(np.floor(spec.wave.data.min()), np.ceil(spec.wave.data.max()), 6)
        wave_interp = np.append(wave_interp, wave_interp[-1] + 6) if wave_interp[-1] < np.ceil(spec.wave.data.max()) else wave_interp

        name_arr = data['FILENAME']
        disp_arr = data['DISPERSER']
        cenwave_arr = data['CENWAVE']
        lp_arr = data['LIFE_ADJ']
        min_arr = data['MINWAVE']
        max_arr = data['MAXWAVE']
        # idcs_no = np.where(data['CENWAVE'] != '3000')[0]
        idcs_no = np.arange(name_arr.size)
        fwhm_cube = np.full((321, wave_interp.size, name_arr.size), np.nan)

        # Treating resolution
        counter = 0
        for i, idx in enumerate(idcs_no):
            disp = disp_arr[idx]
            if disp != 'G185M':
                lsf_path = lsf_folder/f'{obj}_aa_LSFTable_{disp}_{cenwave_arr[idx]}_LP{int(lp_arr[idx])}_cn_convolved.dat'
            else:
                lsf_path = lsf_folder/f'{obj}_nuv_model_lsf_convolved.dat'

            # Open file and generate the interpolator
            matrix_arr = np.loadtxt(lsf_path)
            wave_lsf, lsf_matrix = matrix_arr[0, :], matrix_arr[1:, :]
            #
            # fig, ax = plt.subplots(figsize=(10, 5))
            # pm = ax.pcolormesh(wave_lsf, np.arange(lsf_matrix.shape[0]), np.ma.masked_invalid(lsf_matrix))
            # ax.set_xlabel('Wavelength (Å)', fontsize=18)
            # ax.set_ylabel('LSF', fontsize=18)
            # fig.colorbar(pm, ax=ax)
            # plt.tight_layout()
            # plt.show()

            idcs_x = (wave_interp > wave_lsf.min()) & (wave_interp < wave_lsf.max())
            idcs_y = (110, 211) if disp == 'G185M' else (0, 322)

            f = interp1d(wave_lsf, lsf_matrix, axis=1, kind='linear')  # or 'cubic', 'nearest', etc.
            fwhm_interp = f(wave_interp[idcs_x])

            fwhm_cube[idcs_y[0]:idcs_y[1], idcs_x, i] = fwhm_interp
            counter += 1

        obj_hasp_lsf = np.nanmean(fwhm_cube, axis=2)

        # Set nan values to zero
        obj_hasp_lsf[np.isnan(obj_hasp_lsf)] = 0

        # Remove columns which are zero
        idcs_non0 = obj_hasp_lsf.any(axis=0)
        obj_hasp_lsf = obj_hasp_lsf[:, idcs_non0]
        wave_interp = wave_interp[idcs_non0]

        fig, ax = plt.subplots(figsize=(10,5))
        pm = ax.pcolormesh(wave_interp, np.arange(obj_hasp_lsf.shape[0]), np.ma.masked_invalid(obj_hasp_lsf))
        ax.set_xlabel('Wavelength (Å)', fontsize=18)
        ax.set_ylabel('LSF', fontsize=18)
        ax.set_title(f'{obj} HASP spectra mean LSF ({counter} spectra (CENWAVE {np.unique(cenwave_arr[idcs_no])}))', fontsize=20)
        fig.colorbar(pm, ax=ax)
        plt.tight_layout()
        plt.savefig(lsf_obj_folder/f'{obj}_hasp_lsf.png')
        # plt.show()

        matrix_with_row = np.vstack([wave_interp, obj_hasp_lsf])
        np.savetxt(lsf_obj_folder/f'{obj}_hasp_lsf.txt', matrix_with_row)