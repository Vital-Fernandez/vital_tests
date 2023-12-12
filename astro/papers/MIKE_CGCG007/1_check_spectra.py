import numpy as np
import lime
from pathlib import Path
from tools import mike_load_function

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

# Create sample wth the nights and the orders
file_log = lime.load_log('index_log.txt', levels=['id', 'file'])
sample_mike = lime.Sample(file_log, load_function=mike_load_function, redshift=z_obj)

# Loop through the spectra and bands
for i, spec_path in enumerate(spectra_list):
    print(f'{spec_path}')
    idcs = sample_mike.log.file_path == spec_path.as_posix()
    sample_mike.plot.spectra(idcs, plt_cfg=fig_conf, legend_handle=None)

for order in np.arange(33):

    idcs = (sample_mike.log.order == order) & (sample_mike.log.arm == "blueArm")
    sample_mike.plot.spectra(idcs, plt_cfg=fig_conf, maximize=True)