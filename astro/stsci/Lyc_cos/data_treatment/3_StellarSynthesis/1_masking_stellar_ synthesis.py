import re
import lime
import numpy as np
from pathlib import Path
from astro.stsci.tools import unpack_lines
from matplotlib.widgets import SpanSelector
from matplotlib import pyplot as plt
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
lyman_alpha_folder = project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'
lsf_obj_folder = project_folder/'LyC_leakers_COS'/'lsf_objects'

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')
grating_list = ['G130M', 'G160M', 'G185M']
lime_voigtfit_conv = cfg_sample['LiMe_voigtfit_conversion']

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Order for objects
target_list = sample_df.object.unique()
target_list = ['SBS0335052', 'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
               'UGCA281', 'NGC1705',  'NGC4861',  'Haro2', 'He2-10', 'MRK1450',         # Codo
               'UGC4483', 'VIIZw403',  'NGC2366',                                       # Dent
               'Haro11_A', 'Haro11_B', 'Haro11_C',
               'IZw18_SE']

# target_list = ['Haro2', 'He2-10']

# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    if i >= 0:

        # Object folders the files
        print(f'{i}) Galaxy {obj}, z = {cfg_sample['Galaxy_redshifts'][obj]}')
        input_folder_single = obs_folder / 'LyC_leakers_COS' / 'objects_x1d' / f'{obj}'
        output_folder_single = obs_folder / 'LyC_leakers_COS' / 'obj_hasp' / f'{obj}'
        opacityLymanAlpha_df_path = lyman_alpha_folder/f'{obj}_LyAlpha_lines_frame.txt'
        opacity_df_path = output_folder/f'{obj}_metals_lines_frame.txt'
        voigfitz_reg = output_folder/f'{obj}_metals_best_fit.reg'
        obj_lsf_file = f"{lsf_obj_folder}/{obj}_hasp_lsf.txt"
        metals_cfg = cfg_sample['voigtfit_metals_params'][obj]

        # Object folders the files
        spec = lime.Spectrum.from_file(output_folder/f'{obj}_rebinned_LymanAlphaCorrected.txt', instrument='text')

        # Sort the line selection line lists
        list_science, list_masked = unpack_lines(spec, metals_cfg, science_groups=['ISM'],
                                                 mask_groups=['airglow', 'MW', 'nebular', 'stellar'])

        fig_cfg = {"legend.fontsize": 10}
        spec.plot.spectrum(line_list=list_masked + list_science, in_fig=None, fig_cfg=fig_cfg, ax_cfg={'title': obj})
        span = SpanSelector(spec.plot.ax, onselect, direction='horizontal', useblit=True, button=3,
                            interactive=True, props=dict(alpha=0.3))
        spec.plot.ax.set_xlim(spec.wave.min(), spec.wave.max())
        plt.tight_layout()
        spec.plot.show()