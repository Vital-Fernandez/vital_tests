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


lime.theme.set_style('dark')

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
               'Haro11_A', 'Haro11_B', 'Haro11_C',
               'IZw18_SE']

# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    if i == 17:

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
        # Rebinned
        spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)

        # Divide by the HI profile
        norm_profile = absorption_spectrum(spec, ['H1_1216A'], HI_cfg, plot_fit=False)
        spec.flux = spec.flux /norm_profile
        spec.err_flux = spec.err_flux /norm_profile

        # Normalization
        spec = spec.retrieve.normalization(**metals_cfg['normalization'])

        # Spectra masking
        if 'Metals_retrieve' in metals_cfg:
            spec = spec.retrieve.spectrum(**metals_cfg['Metals_retrieve'])

        # Sort the line selection line lists
        list_science, list_masked = unpack_lines(spec, metals_cfg, science_groups=['ISM'],
                                                 mask_groups=['airglow', 'MW', 'nebular', 'stellar'])

        # Get line list extra lines
        bands_target = spec.retrieve.lines_frame(band_vsigma=70, **metals_cfg['LyC_lines'])
        # final_list = metals_cfg['LyC_lines']['line_list'] + metals_cfg['LyC_lines']['grouped_lines']
        # bands_target = bands_target.loc[bands_target.index.isin(final_list)]

        if 'clean_lines' in metals_cfg:
            bands_target = bands_target[~bands_target.index.isin(metals_cfg['clean_lines'])]

        # bands_target.drop('Si3_1206A', inplace=True)
        # Run Voigtfit
        spec.plot.spectrum(bands=bands_target, line_list=list_masked + list_science, rest_frame=False)
        run_VoigtFit(output_folder/f'{obj}_metals', spec, bands_target, fit_cfg=metals_cfg['LyC_lines']['fit_cfg'],
                     conv_dict=lime_voigtfit_conv, output_toml=output_folder/'metals_results.toml',
                     obj_redshift=spec.redshift, voigt_default_params=cfg_sample['voigt_default_params'])



'''
XXXXX.LyC_lines.map_band_vsigma.S2_1251A = 70
XXXXX.LyC_lines.map_band_vsigma.S2_1251A_o-MW = 70
XXXXX.LyC_lines.map_origin.none.rejected_lines = ['C3_1176A']
XXXXX.LyC_lines.map_origin.MW.rejected_lines = ['C3_1176A']

XXXXX.LyC_lines.map_origin.none.line_list = ['Fe2_1097A',
				             'Fe2_1142A', 'Fe2_1143A','Fe2_1145A',
				             'S3_1190A', 'Si2_1190A', 'Si2_1193A',
				             'Si3_1206A',
				             'S2_1251A', 'S2_1254A',
				             'S2_1260A', 'Fe2_1260A', 'Si2_1260A',
				             'Si2_1304A',
				             'C2_1335A', 'C2_1336A', 'C2_1335.7A',
				             'Si2_1527A',
				             'Fe2_1608A',
				             'Al2_1671A']

XXXXX.LyC_lines.map_origin.MW.line_list = ['Al2_1671A']


XXXXX.LyC_lines.fit_cfg.Fe2_1142A_b = 'Fe2_1142A+Fe2_1143A+Fe2_1145A+Fe2_1142A_o-MW+Fe2_1143A_o-MW+Fe2_1145A_o-MW'
XXXXX.LyC_lines.fit_cfg.S3_1190A_b = 'S3_1190A+Si2_1190A+Si2_1193A+S3_1190A_o-MW+Si2_1190A_o-MW+Si2_1193A_o-MW'
XXXXX.LyC_lines.fit_cfg.S2_1251A_b = 'S2_1251A+S2_1254A'
XXXXX.LyC_lines.fit_cfg.S2_1260A_b = 'S2_1260A+Fe2_1260A+Si2_1260A'
XXXXX.LyC_lines.fit_cfg.C2_1335A_b = 'C2_1335A+C2_1336A+C2_1335.7A'

XXXXX.LyC_lines.grouped_lines = ['Fe2_1142A_b', 'S3_1190A_b', 'Si3_1206A_b', 'S2_1251A_b', 'S2_1260A_b',
                                 'C2_1335A_b', 'Si2_1527A_b', 'Fe2_1608A_b']
'''
