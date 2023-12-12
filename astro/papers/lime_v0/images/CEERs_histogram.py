import numpy as np
import pandas as pd
import lime
from pathlib import Path
from matplotlib import pyplot as plt, rc_context

# Sample file
cfg = lime.load_cfg(Path(r'D:\Pycharm Projects\CEERs_field\reduction_v0.7\CEERs_DR0.7.toml'))
sample_name = cfg['file_structure']['release']
norm_flux = cfg['sample_data']['norm_flux']
fits_folder = Path(cfg['file_structure']['spectra_folder'])
bands_folder = Path(cfg['file_structure']['bands_folder'])/sample_name
data_folder = Path(cfg['file_structure']['data_folder'])
logs_folder = Path(cfg['file_structure']['logs_folder'])

# File log dataframe
file_log_address = data_folder/f'{sample_name}_files_log.txt'
file_sample = lime.Sample(file_log_address, levels=["sample", "id", "file"])

# File log dataframe
flux_log_address = data_folder/f'{sample_name}_flux_log.txt'
flux_df = lime.load_log(flux_log_address, levels=["sample", "id", "file", "line"])


z_df = lime.redshift_calculation(flux_df, weight_parameter='gauss_flux')

# Check for duplicated values in the 'Level1' index and create a boolean mask
duplicated_mask = z_df.index.get_level_values('id').duplicated(keep='first')
z_df_unique = z_df[~duplicated_mask]

n_obs = z_df.index.size
n_unique = z_df.index.get_level_values('id').unique().size
lime.save_log(z_df, data_folder/f'{sample_name}_redshift_log.txt')
print(f'Number of observations: {n_obs}, with {n_unique} objects')

fig_cfg = {'figure.dpi': 300,
           'figure.figsize': (10, 8),
           'axes.titlesize': 14,
           'axes.labelsize': 22,
           'legend.fontsize': 12,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20}



with rc_context(fig_cfg):
    fig, ax = plt.subplots(dpi=300)
    ax.hist(z_df_unique.z_mean.to_numpy(), color='skyblue', edgecolor='black')  # Adjust the number of bins as needed
    plt.xlabel('Redshift')
    plt.ylabel('Unique object count')

    # Show the histogram
    plt.tight_layout()
    # plt.show()
    plt.savefig(f'{sample_name}_redshift_histogram.png')

