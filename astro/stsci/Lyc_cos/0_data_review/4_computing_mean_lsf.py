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
cfg_sample = lime.load_cfg('../Lyc_cos.toml')

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

    if i >= 0:

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
        idcs_no = np.where(data['CENWAVE'] != '3000')[0]
        fwhm_cube = np.full((321, wave_interp.size, name_arr.size), np.nan)

        # Treating resolution
        counter = 0
        for i, idx in enumerate(idcs_no):
            disp = disp_arr[idx]
            print(disp, cenwave_arr[idx], int(lp_arr[idx]))
            if disp != 'G185M':
                lsf_path = lsf_folder/f'{obj}_aa_LSFTable_{disp}_{cenwave_arr[idx]}_LP{int(lp_arr[idx])}_cn_convolved.dat'
            else:
                lsf_path = lsf_folder/f'{obj}_nuv_model_lsf_convolved.dat'

            # Open file and generate the interpolator
            matrix_arr = np.loadtxt(lsf_path)
            wave_lsf, lsf_matrix = matrix_arr[0, :], matrix_arr[1:, :]

            idcs_x = (wave_interp > wave_lsf.min()) & (wave_interp < wave_lsf.max())
            idcs_y = (110, 211) if disp == 'G185M' else (0, 322)

            f = interp1d(wave_lsf, lsf_matrix, axis=1, kind='linear')  # or 'cubic', 'nearest', etc.
            fwhm_interp = f(wave_interp[idcs_x])

            fwhm_cube[idcs_y[0]:idcs_y[1], idcs_x, i] = fwhm_interp
            counter += 1

        obj_hasp_lsf = np.nanmean(fwhm_cube, axis=2)

        obj_hasp_lsf[np.isnan(obj_hasp_lsf)] = 0

        # spec.plot.spectrum()
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

        # for j in np.arange(obj_hasp_lsf.shape[1]):
        #     fig, ax = plt.subplots()
        #     ax.plot(obj_hasp_lsf[:, j], label=f'{i}')
        #     ax.legend()
        #     plt.show()

        # for i in np.arange(fwhm_cube.shape[2]):
        #     for j in np.arange(fwhm_cube.shape[1]):
        #         fig, ax = plt.subplots()
        #         ax.plot(fwhm_cube[:, j, i], label=f'{i}, {j}')
        #         ax.legend()
        #         plt.show()
# # Wavelength range for interpolation
        # min_G130M = np.min(min_arr[np.where(disp_arr == 'G130M')[0]])
        #
        #
        # if 'G160M' in disp_arr:
        #     max_G160M = np.max(max_arr[np.where(disp_arr == 'G160M')[0]])
        # else:
        #     max_G160M = np.max(min_arr[np.where(disp_arr == 'G130M')[0]])
        #
        # if 'G185M' in disp_arr:
        #     min_G180M = np.min(min_arr[np.where(disp_arr == 'G185M')[0]])
        #     # min_G180M = min_G180M if min_G180M > max_G160M else max_G160M
        #     max_G180M = np.max(max_arr[np.where(disp_arr == 'G185M')[0]])
        # else:
        #     min_G180M, max_G180M = None, None
        #
        # print(obj, min_G130M, max_G160M, min_G180M, max_G180M)

        # for m, disp in enumerate(np.unique(disp_arr)):
        #     idcs_disp = np.where(disp_arr == disp)[0]
        #
        #     for idx in idcs_disp:
        #         if disp != 'G185M':
        #             lsf_path = lsf_folder/f'{obj}_aa_LSFTable_{disp}_{cenwave_arr[idx]}_LP{int(lp_arr[idx])}_cn_convolved.dat'
        #         else:
        #             lsf_path = lsf_folder/f'{obj}_nuv_model_lsf_convolved.dat'
        #         matrix_arr = np.loadtxt(lsf_path)
        #         wave_lsf, lsf_matrix = matrix_arr[0, :], matrix_arr[1:, :]
        #         # print(disp, cenwave_arr[idx], lp_arr[idx], f'{min_arr[idx]:0.2f}', f'{max_arr[idx]:0.2f}', np.mean(np.diff(wave_lsf)))
        #         print(disp, cenwave_arr[idx], lp_arr[idx], f'{wave_lsf[0]}', f'{wave_lsf[-1]}', np.mean(np.diff(wave_lsf)))


