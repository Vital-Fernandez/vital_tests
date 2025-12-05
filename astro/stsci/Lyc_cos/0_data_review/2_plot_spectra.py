import re
from pathlib import Path
import numpy as np
import lime
from astropy.io import fits
from matplotlib import pyplot as plt

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))

idcs_manual = sample_df.index.get_level_values('state') == 'aspec_manual'
obj_list = sample_df.loc[idcs_manual, 'object'].unique()

for i, obj_name in enumerate(obj_list):

    # Get subspec
    idcs_cspec = (sample_df.object == obj_name) & (sample_df.index.get_level_values('state') == 'cspec_manual')
    sub_sample = sample_df.loc[idcs_cspec]

    # Manual cspec
    if sub_sample.index.size == 0:
        print(f'-{obj_name} NO SPECTRA')

    else:
        redshift = cfg_sample['Galaxy_redshifts'][obj_name]
        file_path = obs_folder / sample_df.loc[idcs_manual & (sample_df.object == obj_name), 'filepath'].values[0]
        print(f'- {obj_name} at z = {redshift} ({file_path.stem})')

        spec = lime.Spectrum.from_file(file_path, instrument='COS', redshift=redshift)
        spec.plot.spectrum(in_fig=None, label=f'{file_path.stem}')

        cmap = plt.get_cmap("plasma_r")
        color_list = cmap(np.linspace(0, 1, len(sub_sample)))

        for j, idx in enumerate(sub_sample.index):
            file_path = obs_folder / sub_sample.loc[idx, 'filepath']
            c_spec = lime.Spectrum.from_file(file_path, instrument='COS', redshift=redshift)
            spec.plot.ax.step(c_spec.wave, c_spec.flux, where='mid', label=f'{file_path.stem}', color=color_list[j], linestyle='--')

        # Plot wording
        spec.plot.ax.legend()
        spec.plot.ax.set_title(obj_name)
        plt.tight_layout()
        plt.show()

        # for j, idx in enumerate(sub_sample.index):
        #     file_path = obs_folder/sub_sample.loc[idx, 'filepath']
        #
        #     if j == 0:
        #         print(fits.info(file_path))
        #         spec = lime.Spectrum.from_file(file_path, instrument='cos', redshift=redshift, norm_flux=1e-17)
        #         label = f'{obj_name} - {idx[2]} - {idx[3]}'
        #         spec.plot.spectrum(label=label, maximize=True, in_fig=None)
        #     else:
        #         spec_j = lime.Spectrum.from_file(file_path, instrument='cos', redshift=redshift, norm_flux=1e-17)
        #         label = f'{obj_name} - {idx[2]} - {idx[3]}'
        #         spec.plot.ax.step(spec_j.wave, spec_j.flux, where='mid', linestyle=':', label=label)
        #     print('--', file_path, idx)
        #
        # # MAST cspec
        # idcs = (sample_df.object == obj_name) & (sample_df.index.get_level_values('state') == 'hasp_mast')
        # sub_sample = sample_df.loc[idcs]
        # n_cspec_mast = sub_sample.index.size
        # if sub_sample.index.size == 0:
        #     print(f'-{obj_name} NO SPECTRA')
        # else:
        #     for j, idx in enumerate(sub_sample.index):
        #         file_path = obs_folder / sub_sample.loc[idx, 'filepath']
        #         spec_j = lime.Spectrum.from_file(file_path, instrument='cos', redshift=redshift, norm_flux=1e-17)
        #         label = f'{obj_name} - {idx[2]} - {idx[3]}'
        #         spec.plot.ax.step(spec_j.wave, spec_j.flux, where='mid', linestyle='--', label=label)
        #
        # spec.plot.ax.legend()
        # spec.plot.ax.set_title(f'{obj_name}, n_x1d = {n_objects}, n_cspec_manual={n_cspec_manual} , n_cspec_mast = {n_cspec_mast}')
        # spec.plot.ax.set_yscale('symlog', linthresh=10)
        # plt.tight_layout()
        # plt.show()