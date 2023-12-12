import numpy as np
import pandas as pd

import lime
from pathlib import Path
from tools import read_mike_fits

cfg = lime.load_cfg('mike_cgcg007.toml')
data_folder = Path(cfg['folders']['data'])
treatment_folder = Path(cfg['folders']['treatment'])
parent_bands = Path(cfg['folders']['parent_bands_file'])
redshift_path = Path(cfg['folders']['redshift_order_file'])

file_df = pd.DataFrame(columns=['id', 'file', 'night', 'arm', 'order', 'wmin', 'wmax', 'file_path'])
file_df.set_index(['id', 'file'], inplace=True)

# Spectra list
spectra_list = list(data_folder.glob('*.fits'))

# Observation dat
z_obj = cfg['CGCG007_025']['redshift']

# Plot configuration
fig_conf = {'figure.dpi': 200, 'lines.linewidth': 0.75, 'font.size': 8}

# Loop through the spectra and bands
for i, spec_path in enumerate(spectra_list):

    print(f'{spec_path}')
    wave_m, flux_m, err_m = read_mike_fits(spec_path)

    for j_order in np.arange(wave_m.shape[0]):

        obs_ref = f'{spec_path.stem}_{j_order}'
        print(f'- {obs_ref}')

        wave_array, flux_array, err_array = wave_m[j_order, :], flux_m[j_order, :], err_m[j_order, :]

        id_label = f'{spec_path.stem}_{j_order}'
        file = spec_path.name
        comps = spec_path.name.split('_')

        file_df.loc[(id_label, file), 'night'] = int(comps[2][5:-5])
        file_df.loc[(id_label, file), 'arm'] = comps[1]
        file_df.loc[(id_label, file), 'order'] = j_order
        file_df.loc[(id_label, file), 'file_path'] = spec_path.as_posix()

        file_df.loc[(id_label, file), 'wmin'] = wave_array[0]
        file_df.loc[(id_label, file), 'wmax'] = wave_array[-1]
        # spec_j = lime.Spectrum(wave_array, flux_array, err_array, redshift=z_obj)
        # spec_j.plot.spectrum(rest_frame=False, fig_cfg=fig_conf)

        # Check redshift
        # spec_j.check.redshift(obs_ref, parent_bands, redshift_path)

        # Review bands
        # order_bands_file = treatment_folder/f'{obs_ref}_bands.txt'
        # spec_j.check.bands(order_bands_file, parent_bands, object_ref=obs_ref, redshift_log=redshift_path,
        #                    maximize=True)

lime.save_log(file_df, f'index_log.txt')