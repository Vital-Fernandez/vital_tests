import re
import lime
from pathlib import Path
import numpy as np
from astro.stsci.tools import run_VoigtFit

lime.theme.set_style('dark')

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
output_folder = project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'

# Cfg file
cfg_sample = lime.load_cfg('../Lyc_cos.toml')
grating_list = ['G130M', 'G160M', 'G185M']
lime_voigtfit_conv = cfg_sample['LiMe_voigtfit_conversion']

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Order for objects
target_list = sample_df.object.unique()
target_list = ['SBS0335052', 'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
               'UGCA281', 'NGC1705',  'NGC4861',  'Haro2', 'He2-10', 'MRK1450',
               'UGC4483', 'VIIZw403',  'NGC2366',
               'Haro11_A', 'Haro11_B', 'Haro11_C',
               'IZw18_SE']

# Loop through the targets and generate the HASP files
check_masks = False
for i, obj in enumerate(target_list):

    if i >= 0:

        # Object folders the files
        print(f'{i}) Galaxy {obj}, z = {cfg_sample['Galaxy_redshifts'][obj]}')
        input_folder_single = obs_folder / 'LyC_leakers_COS' / 'objects_x1d' / f'{obj}'
        output_folder_single = obs_folder / 'LyC_leakers_COS' / 'obj_hasp' / f'{obj}'
        LyA_cfg = cfg_sample['voigtfit_lyAlpha_params'][obj]

        # Index the files
        idcs_out = (sample_df.object == obj) & (sample_df.index.get_level_values('state') == 'aspec_manual')

        # Load the files
        pname = obs_folder / sample_df.loc[idcs_out].filepath[0]
        spec = lime.Spectrum.from_file(pname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj],
                                       **LyA_cfg['from_file'])

        # Rebinned
        spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)
        spec.fit.continuum()

        # Normalization
        spec = spec.retrieve.normalization(**LyA_cfg['normalization'])

        # Spectra masking
        spec = spec.retrieve.spectrum(**LyA_cfg['LyA_retrieve'])
        spec.save_spectrum(output_folder/f'{obj}_LyAlpha_spec.txt')

        # Get the lines
        line_df = spec.retrieve.lines_frame(fit_cfg=LyA_cfg['fit_cfg'], **LyA_cfg['LyA_lines_frame'])

        # Run Voigtfit
        # spec.plot.spectrum(bands=line_df, ax_cfg={'title': f'{obj}'})
        run_VoigtFit(output_folder/f'{obj}_LyAlpha', spec, line_df, fit_cfg=LyA_cfg['fit_cfg'],
                     conv_dict=lime_voigtfit_conv, output_toml=output_folder/'LyAlpha_results.toml',
                     obj_redshift=spec.redshift, voigt_default_params=cfg_sample['voigt_default_params'], var_z=False)



