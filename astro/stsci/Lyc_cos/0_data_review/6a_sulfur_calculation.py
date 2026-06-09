import re
import lime
from pathlib import Path
from astro.stsci.tools import run_VoigtFit, unpack_lines, absorption_spectrum, add_opacity_profile


lime.theme.set_style('dark')

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
output_folder = project_folder/'LyC_leakers_COS'/'metals_line_fitting'
lyman_alpha_folder = project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'
lsf_obj_folder = project_folder/'LyC_leakers_COS'/'lsf_objects'

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

    if i == 15:

        # Object folders the files
        print(f'{i}) Galaxy {obj}, z = {cfg_sample['Galaxy_redshifts'][obj]}')
        input_folder_single = obs_folder / 'LyC_leakers_COS' / 'objects_x1d' / f'{obj}'
        output_folder_single = obs_folder / 'LyC_leakers_COS' / 'obj_hasp' / f'{obj}'
        opacityLymanAlpha_df_path = lyman_alpha_folder/f'{obj}_LyAlpha_lines_frame.txt'
        opacity_df_path = output_folder/f'{obj}_metals_lines_frame.txt'
        obj_lsf_file = f"{lsf_obj_folder}/{obj}_hasp_lsf.txt"

        HI_cfg = lyAlpha_data[f'{obj}_LyAlpha_results']
        metals_cfg = cfg_sample['voigtfit_metals_params'][obj]

        if obj in cfg_sample['voigtfit_sulfur_params']:
            sulfur_cfg = cfg_sample['voigtfit_sulfur_params'][obj]
            voigfit_reg = output_folder / f'{obj}_sulfur_best_fit.reg'
            opacity_df_path = output_folder / f'{obj}_sulfur_lines_frame.txt'

            # Index the files
            idcs_out = (sample_df.object == obj) & (sample_df.index.get_level_values('state') == 'aspec_manual')

            fname = '/home/vital/Astrodata/STScI/LyC_leakers_COS/obj_hasp/Haro11_A/hst_haro11-a_aspec.fits'
            spec = lime.Spectrum.from_file(fname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj])
            spec = lime.Spectrum.from_file(fname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj])

            # Load the files
            pname = obs_folder / sample_df.loc[idcs_out].filepath[0]
            spec = lime.Spectrum.from_file(pname, instrument='cos', redshift=cfg_sample['Galaxy_redshifts'][obj],
                                           **metals_cfg['from_file'])
            # Rebinned
            spec = spec.retrieve.rebinned(pixel_number=6, return_spectrum=True)

            # Divide by the HI profile
            norm_Lyman_alpha = absorption_spectrum(opacityLymanAlpha_df_path, spec)
            norm_Lyman_alpha[norm_Lyman_alpha < 0.01] = 1
            spec.flux = spec.flux /norm_Lyman_alpha
            spec.err_flux = spec.err_flux /norm_Lyman_alpha

            # Normalization
            spec = spec.retrieve.normalization(**metals_cfg['normalization'])

            # Spectra masking
            if 'Metals_retrieve' in metals_cfg:
                spec = spec.retrieve.spectrum(**metals_cfg['Metals_retrieve'])

            # Save the spectrum
            spec.save_spectrum(fname=output_folder/f'{obj}_metals_spec.txt')

            # Sort the line selection line lists
            list_science, list_masked = unpack_lines(spec, metals_cfg, science_groups=['ISM'],
                                                     mask_groups=['airglow', 'MW', 'nebular', 'stellar'])

            # Get line list extra lines
            bands_target = spec.retrieve.lines_frame(band_vsigma=70, **sulfur_cfg['LyC_lines'])

            # Plot inputs
            fig_cfg = {"legend.fontsize": 10,}
            spec.plot.spectrum(bands=bands_target, line_list=list_masked + list_science, rest_frame=False, in_fig=None, fig_cfg=fig_cfg, ax_cfg={'title': obj})
            add_opacity_profile(spec, opacity_df_path, voigtfit_pname=voigfit_reg)

            # Run Voigtfit
            run_VoigtFit(output_folder/f'{obj}_sulfur', spec, bands_target, fit_cfg=sulfur_cfg['LyC_lines']['fit_cfg'],
                         conv_dict=lime_voigtfit_conv, output_toml=output_folder/'metals_results.toml',
                         obj_redshift=spec.redshift, voigt_default_params=cfg_sample['voigt_default_params'],
                         lsf_file=obj_lsf_file)
