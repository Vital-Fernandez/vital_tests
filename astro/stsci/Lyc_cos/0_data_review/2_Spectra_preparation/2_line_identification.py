import re
import lime
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np

def on_click(event):
    toolbar = event.canvas.toolbar

    # If a tool (pan/zoom) is active, do nothing
    if toolbar is not None and toolbar.mode != '':
        return

    if event.button == 3 and event.inaxes is not None:
        print(f"x = {event.xdata:.3f}")

lime.theme.set_style('dark')

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')
grating_list = ['G130M', 'G160M', 'G185M']

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Order for objects
target_list = sample_df.object.unique()
target_list = ['SBS0335052',  'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
               'UGCA281', 'NGC1705',  'NGC4861',  'Haro2', 'He2-10', 'MRK1450',         # Codo
               'UGC4483', 'VIIZw403',  'NGC2366',                                       # Dent
               'Haro11_A', 'Haro11_B', 'Haro11_C',
               'IZw18_SE']


# N I, N II, C II, C III, Fe II, Fe III, O I, O II, Si II, Si III, Si IV, S II, S III, Al II, N IV, N V
# Missing line 1137.467304
# Missing line 1474.467304

# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):
    if i >= 0:

        # Object folders the files
        print(f'{i}) Galaxy {obj}, z = {cfg_sample['Galaxy_redshifts'][obj]}')
        input_folder_single = obs_folder / 'LyC_leakers_COS' / 'objects_x1d' / f'{obj}'
        output_folder_single = obs_folder / 'LyC_leakers_COS' / 'obj_hasp' / f'{obj}'

        # Index the files
        idcs_target = (sample_df.object == obj)
        idcs_in = idcs_target & (sample_df.index.get_level_values('state') == 'x1d') & (sample_df.grating.isin(grating_list))
        idcs_out = idcs_target & (sample_df.index.get_level_values('state') == 'aspec_manual')

        # Load the files
        pname = obs_folder / sample_df.loc[idcs_out].filepath[0]
        spec = lime.Spectrum.from_file(pname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj], crop_flux=(0, 99))

        # Rebinned
        spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)


        spec = lime.Spectrum.from_file(file_path, instrument='COS', redshift=redshift)
        spec_rebinn = spec.retrieve.rebinned(pixel_number=6, constant_pixel_width=True, return_spectrum=True)
        spec_rebinn.plot.spectrum()

        # # Create horizontal SpanSelector
        # ax_cfg = {'title': f'Galaxy {obj}'}
        # fig_cfg = {"legend.fontsize":  10}
        # spec.plot.spectrum(in_fig=None, lines=['ISM', 'Nebular', 'Airglow', "Milky_way"], rest_frame=True,
        #                    fit_cfg=cfg_sample, obj_cfg_prefix=obj,
        #                    fig_cfg=fig_cfg, ax_cfg=ax_cfg)
        # spec.plot.fig.canvas.mpl_connect("button_press_event", on_click)
        # plt.tight_layout()
        # spec.plot.show()

# 1142.36565
# 1143.22573
# 1144.93920