import re
import lime
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from lime.fitting.lines import c_KMpS



lime.theme.set_style('dark')

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
bands_folder = project_folder/'LyC_leakers_COS/bands'
fwhm_folder = project_folder/'LyC_leakers_COS/acq_im_fwhm'
lines_folder = project_folder/'LyC_leakers_COS/line_measurements'

# Cfg file
cfg_sample = lime.load_cfg('../Lyc_cos.toml')
grating_list = ['G130M', 'G160M', 'G185M']

# Sample file
sample_df = lime.load_frame(project_folder / 'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Order for objects
target_list = sample_df.object.unique()
target_list = ['SBS0335052', 'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
               'UGCA281', 'NGC1705', 'NGC4861', 'Haro2', 'He2-10', 'MRK1450',  # Codo
               'UGC4483', 'VIIZw403', 'NGC2366',  # Dent
               'Haro11_A', 'Haro11_B', 'Haro11_C',
               'IZw18_SE']

# Loop through the targets and generate the HASP files
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
        spec = lime.Spectrum.from_file(pname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj],
                                       crop_flux=(0,99))

        # Rebinned
        spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)


        # Define the line bands for the Galaxy
        lines_neb = cfg_sample[f'default_line_fitting']['Nebular']['labels']
        bands_nebular = spec.retrieve.lines_frame(line_list=lines_neb)
        file_name = bands_folder/f'{obj}_nebular_bands.txt'

        # Manually review the bands
        # spec.check.bands(file_name, line_list=lines_neb, show_continua=True)

        # Perform the fittings
        spec.fit.frame(file_name, cfg_sample, obj_cfg_prefix=obj, cont_source='adjacent')
        spec.plot.grid(lines_folder/f'{obj}_nebular_lines_grid.png')
        spec.save_frame(lines_folder/f'{obj}_nebular_lines_log.txt')

        df = spec.frame
        w3, w4 = df['w3'] * (1 + spec.redshift), df['w4'] * (1 + spec.redshift)
        idx3, idx4 = np.searchsorted(spec.wave_rest, (w3, w4))
        res_arr = (w4 - w3) / (idx4- idx3)
        FWHM_pixels = df['FWHM_p']/res_arr
        lime.save_frame(fwhm_folder/f'{obj}_lines_FWHM_pixels.txt', FWHM_pixels)

        # # Create horizontal SpanSelector
        ax_cfg = {'title': f'Galaxy {obj}'}
        fig_cfg = {"legend.fontsize": 10}
        spec.plot.spectrum(bands=bands_nebular, lines=['ISM', 'Nebular', 'Airglow', "Milky_way"],
                           fit_cfg=cfg_sample, obj_cfg_prefix=obj, fig_cfg=fig_cfg, ax_cfg=ax_cfg, rest_frame=False,
                           show_err=True)

