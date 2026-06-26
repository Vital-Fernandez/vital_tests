import re

import numpy as np
from Cython.Shadow import returns

import lime
from pathlib import Path
from astro.stsci.tools import list_files_with_extension, add_cos_obs, run_hasp_wrapper, move_files
import lineid_plot
from matplotlib.widgets import SpanSelector
import matplotlib.pyplot as plt
lime.theme.set_style('dark')

Table_lists = []
def onselect(xmin, xmax):
    Table_lists.append(np.round([xmin, xmax], 2).tolist())
    print(f"\nSelections: {Table_lists}")
    # print(f"Selected x-range: ({xmin:.2f}, {xmax:.2f}),")


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Get list of targets
target_list = sample_df.object.unique()
grating_list = ['G130M', 'G160M', 'G185M']

# Lines list
line_labels = cfg_sample['line_properties']['emis_lines']
particles, trans_arr, latex_arr = lime.label_decomposition(line_labels)

target_list = ['SBS0335052',  'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
               'UGCA281', 'NGC1705',  'NGC4861',  'Haro2', 'He2-10', 'MRK1450',# Codo
               'UGC4483', 'VIIZw403',  'NGC2366',          # Dent
               'Haro11_A', 'Haro11_B', 'Haro11_C',
               'IZw18_SE']

# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    if i >= 0:

        # Object folders the files
        print(f'{i}) Galaxy {obj}')
        input_folder_single = obs_folder/'LyC_leakers_COS'/'objects_x1d'/f'{obj}'
        output_folder_single = obs_folder/'LyC_leakers_COS'/'obj_hasp'/f'{obj}'

        # Plot the data
        idcs_target = (sample_df.object == obj)
        idcs_in = idcs_target & (sample_df.index.get_level_values('state') == 'x1d') & (sample_df.grating.isin(grating_list))
        idcs_out = idcs_target & (sample_df.index.get_level_values('state') == 'aspec_manual')

        # Load the spectrum
        pname = obs_folder / sample_df.loc[idcs_out].filepath[0]
        spec = lime.Spectrum.from_file(pname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj])
        spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)

        # # Create horizontal SpanSelector
        # spec.plot.spectrum(in_fig=None, rest_frame=False)
        # span = SpanSelector(spec.plot.ax, onselect, direction='horizontal', useblit=True, button=3, interactive=True,
        #                     props=dict(alpha=0.3))
        # plt.tight_layout()
        # spec.plot.show()

        # spec.fit.continuum(degree_list=[4, 5, 6, 6],
        #                    emis_threshold=[6, 5, 4, 3],
        #                    exclude_intvls=cfg_sample['continua_mask'][obj],
        #                    smooth_scale=10,
        #                    plot_steps='last',
        #                    rest_intvls=False,
        #                    fig_cfg={"axes.titlesize": 16},
        #                    ax_cfg={'title': f'Galaxy: {obj}'})
        # spec.plot.spectrum(show_cont=True, rest_frame=True)

        spec = spec.retrieve.normalization(degree_list=[3, 3, 3, 3],
                                           emis_threshold=[6, 5, 4, 3],
                                           exclude_intvls=cfg_sample['continua_mask'][obj],
                                           smooth_scale=10,
                                           return_spectrum=True)
        spec.plot.spectrum(show_cont=True, rest_frame=True)

# # Multi PID objects
# obj_list_multID = group_PID[group_PID > 1].index.tolist()
# for i, obj_name in enumerate(obj_list_multID):
#     idcs_x1d = (sample_df.index.get_level_values('state') == 'x1d') & (sample_df.object == obj_name)
#     file_arr = sample_df.loc[idcs_x1d, 'filepath'].to_numpy()
#     move_files(file_arr, obs_folder, obj_folder_list[i])

# # Run wrapper for multiple PI
# if output_folder_mult.exists():
#     shutil.rmtree(output_folder_mult)
# output_folder_mult.mkdir(parents=True, exist_ok=True)
# for obj_folder in obj_folder_list:
#     run_hasp_wrapper(obj_folder, output_folder_mult, cross_program=True)

# list_hasp = list_files_with_extension(output_folder_mult, '_cspec.fits')
# add_cos_obs(list_hasp, sample_df, 'LyC_cos', 'hasp_mult', cfg_sample)
# lime.save_frame(project_folder/'stsci_samples_v1.csv', sample_df)

# list_hasp = list_files_with_extension(output_folder_mult, '_aspec.fits')
# add_cos_obs(list_hasp, sample_df, 'LyC_cos', 'hasp_aspec', cfg_sample)
# lime.save_frame(project_folder/'stsci_samples_v1.csv', sample_df)