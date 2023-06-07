import numpy as np
import lime
from pathlib import Path
from tools import read_mike_fits

cfg = lime.load_cfg('mike_cgcg007.toml')
data_folder = Path(cfg['folders']['data'])
treatment_folder = Path(cfg['folders']['treatment'])
parent_bands = Path(cfg['folders']['parent_bands_file'])
redshift_path = Path(cfg['folders']['redshift_order_file'])

# Spectra list
spectra_list = list(data_folder.glob('*.fits'))

# Observation dat
z_obj = cfg['CGCG007_025']['redshift']

# Plot configuration
fig_conf = {'figure.dpi': 200, 'lines.linewidth': 0.75, 'font.size': 8}

# Loop through the spectra and bands
for i, spec_path in enumerate(spectra_list):

    wave_m, flux_m, err_m = read_mike_fits(spec_path)

    for j_order in np.arange(wave_m.shape[0]):

        obs_ref = f'{spec_path.stem}_{j_order}'
        print(f'- {obs_ref}')

        wave_array, flux_array, err_array = wave_m[j_order, :], flux_m[j_order, :], err_m[j_order, :]
        spec_j = lime.Spectrum(wave_array, flux_array, err_array, redshift=z_obj)
        spec_j.plot.spectrum(rest_frame=False, fig_cfg=fig_conf)

        # Check redshift
        # spec_j.check.redshift(obs_ref, parent_bands, redshift_path)

        # Review bands
        # order_bands_file = treatment_folder/f'{obs_ref}_bands.txt'
        # spec_j.check.bands(order_bands_file, parent_bands, object_ref=obs_ref, redshift_log=redshift_path,
        #                    maximize=True)
