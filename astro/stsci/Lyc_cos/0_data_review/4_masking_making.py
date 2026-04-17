import re
import lime
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
from lime.fitting.lines import c_KMpS
from lime.io import check_fit_conf
from astropy.io import fits
from astro.stsci.tools import (on_click, IntervalSelector, run_VoigtFit, lines_frame_to_lime_lines,
                               unpack_lines, optical_depth_profile, vfdb_df, profile_normflux, absorption_spectrum)
from matplotlib.widgets import SpanSelector

lime.theme.set_style('dark')

Table_lists = []
def onselect(xmin, xmax):
    Table_lists.append(np.round([xmin, xmax], 2).tolist())
    print(f"\nSelections: {Table_lists}")
    # print(f"Selected x-range: ({xmin:.2f}, {xmax:.2f}),")


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
output_folder = project_folder/'LyC_leakers_COS'/'metals_line_fitting'

# Cfg file
cfg_sample = lime.load_cfg('../Lyc_cos.toml')
grating_list = ['G130M', 'G160M', 'G185M']
lime_voigtfit_conv = cfg_sample['LiMe_voigtfit_conversion']
lyAlpha_data = lime.load_cfg(project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'/'LyAlpha_results.toml')
LyC_list = list(cfg_sample['LyC_target_lines'].keys())

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Order for objects
target_list = sample_df.object.unique()
target_list = ['SBS0335052', 'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
               'UGCA281', 'NGC1705',  'NGC4861',  'Haro2', 'He2-10', 'MRK1450',         # Codo
               'UGC4483', 'VIIZw403',  'NGC2366',                                       # Dent
               'Haro11_A', 'Haro11_B', 'Haro11_C']
               # 'IZw18_SE']

# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    if i == 10:

        # Object folders the files
        print(f'{i}) Galaxy {obj}, z = {cfg_sample['Galaxy_redshifts'][obj]}')
        input_folder_single = obs_folder / 'LyC_leakers_COS' / 'objects_x1d' / f'{obj}'
        output_folder_single = obs_folder / 'LyC_leakers_COS' / 'obj_hasp' / f'{obj}'

        HI_cfg = lyAlpha_data[f'{obj}_LyAlpha_results']
        metals_cfg = cfg_sample['voigtfit_metals_params'][obj]

        # Index the files
        idcs_out = (sample_df.object == obj) & (sample_df.index.get_level_values('state') == 'aspec_manual')

        # Load the files
        pname = obs_folder / sample_df.loc[idcs_out].filepath[0]
        spec = lime.Spectrum.from_file(pname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj],
                                       **metals_cfg['from_file'])
        # spec.plot.spectrum()
        # Rebinned
        spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)

        # Divide by the HI profile
        norm_profile = absorption_spectrum(spec, ['H1_1216A'], HI_cfg, plot_fit=False)
        spec.flux = spec.flux /norm_profile
        spec.err_flux = spec.err_flux /norm_profile

        spec.plot.spectrum(in_fig=None, show_cont=True)
        span = SpanSelector(spec.plot.ax, onselect, direction='horizontal', useblit=True, button=3,
                            interactive=True, props=dict(alpha=0.3))
        # spec.plot.ax.set_xlim(spec.wave.min(), spec.wave.max())
        # spec.plot.ax.plot((spec.wave.min(), spec.wave.max()), (1, 1), linestyle='--', color='yellow')
        plt.tight_layout()
        spec.plot.show()

        # Normalization
        spec = spec.retrieve.normalization(**metals_cfg['normalization'])
        spec.plot.spectrum()




        # # Spectra masking
        # if 'Metals_retrieve' in metals_cfg:
        #     spec = spec.retrieve.spectrum(**metals_cfg['Metals_retrieve'])
        #
        # # Sort the line selection line lists
        # list_science, list_masked = unpack_lines(spec, metals_cfg, science_groups=['ISM'],
        #                                          mask_groups=['airglow', 'MW', 'nebular', 'stellar'])
        #
        # bands_target = spec.retrieve.lines_frame(**metals_cfg['LyC_lines'])
        #
        # # Remove extra lines
        # final_list = metals_cfg['LyC_lines']['line_list'] + metals_cfg['LyC_lines']['grouped_lines']
        # bands_target = bands_target.loc[bands_target.index.isin(final_list)]
        # # bands_target.drop('Si3_1206A', inplace=True)
        #
        # # Run Voigtfit
        # spec.plot.spectrum(bands=bands_target, line_list=list_masked + list_science, rest_frame=False)
        # run_VoigtFit(output_folder/f'{obj}_metals', spec, bands_target, fit_cfg=metals_cfg['LyC_lines']['fit_cfg'],
        #              conv_dict=lime_voigtfit_conv, output_toml=output_folder/'metals_results.toml',
        #              obj_redshift=spec.redshift, voigt_default_params=cfg_sample['voigt_default_params'])




# import re
# import numpy as np
# import lime
# from pathlib import Path
# from matplotlib import pyplot as plt
# from matplotlib.widgets import SpanSelector
# from lime.fitting.lines import c_KMpS
# from astro.stsci.tools import on_click, IntervalSelector
#
# lime.theme.set_style('dark')
#
# Table_lists = []
# def onselect(xmin, xmax):
#     Table_lists.append(np.round([xmin, xmax], 2).tolist())
#     print(f"\nSelections: {Table_lists}")
#     # print(f"Selected x-range: ({xmin:.2f}, {xmax:.2f}),")
#
#
# # Data location
# obs_folder = Path('/home/vital/Astrodata/STScI')
# project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
# output_folder = project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'
#
# # Cfg file
# cfg_sample = lime.load_cfg('../Lyc_cos.toml')
# grating_list = ['G130M', 'G160M', 'G185M']
# lime_voigtfit_conv = cfg_sample['LiMe_voigtfit_conversion']
#
# # Sample file
# sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
# pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
# sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]
#
# # Order for objects
# target_list = sample_df.object.unique()
# target_list = cfg_sample['galaxy_sample']['sn_sort_list']
#
# # Loop through the targets and generate the HASP files
# check_masks = False
# for i, obj in enumerate(target_list):
#
#     if i >= 0:
#
#         # Object folders the files
#         print(f'{i}) Galaxy {obj}, z = {cfg_sample['Galaxy_redshifts'][obj]}')
#         input_folder_single = obs_folder / 'LyC_leakers_COS' / 'objects_x1d' / f'{obj}'
#         output_folder_single = obs_folder / 'LyC_leakers_COS' / 'obj_hasp' / f'{obj}'
#         # LyA_cfg = cfg_sample['voigtfit_lyAlpha_params'][obj]
#
#         # Index the files
#         idcs_out = (sample_df.object == obj) & (sample_df.index.get_level_values('state') == 'aspec_manual')
#
#         # Load the files
#         pname = obs_folder / sample_df.loc[idcs_out].filepath[0]
#         spec = lime.Spectrum.from_file(pname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj], id_label=obj)
#                                        # **LyA_cfg['from_file'])
#
#         # Rebinned
#         spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)
#
#         # # Get the lines
#         # line_df = spec.retrieve.lines_frame(**LyA_cfg['LyA_lines_frame'])
#         #
#         # # Normalization
#         # spec = spec.retrieve.normalization(**LyA_cfg['normalization'])
#         #
#         # # Spectra masking
#         # spec = spec.retrieve.spectrum(**LyA_cfg['LyA_retrieve'])
#         #
#         Table_lists = []
#         spec.plot.spectrum(bands=line_df, in_fig=None, show_cont=True)
#         span = SpanSelector(spec.plot.ax, onselect, direction='horizontal', useblit=True, button=3,
#                             interactive=True, props=dict(alpha=0.3))
#         spec.plot.ax.set_xlim(spec.wave.min(), spec.wave.max())
#         spec.plot.ax.plot((spec.wave.min(), spec.wave.max()), (1, 1), linestyle='--', color='yellow')
#         plt.tight_layout()
#         spec.plot.show()
#
#         # # Spectra masking
#         # spec = spec.retrieve.spectrum(**LyA_cfg['LyA_retrieve'])
#         #
#         # # Get the lines
#         # line_df = spec.retrieve.lines_frame(**LyA_cfg['LyA_lines_frame'])
#         # line_list = lines_frame_to_lime_lines(line_df, LyA_cfg['LyA_lines_frame'].get('fit_cfg'))
#         #
#         # # Run Voigtfit
#         # spec.plot.spectrum(bands=line_df)
