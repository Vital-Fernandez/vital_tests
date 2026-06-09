import re
import lime
import pandas as pd
from pathlib import Path
from astro.stsci.tools import unpack_lines, absorption_spectrum, add_opacity_profile


# def join_dict_dataframes(dfs: dict) -> pd.DataFrame:
#     """
#     Concatenate DataFrames from a dictionary, adding the dict key as an
#     additional first level in a MultiIndex.
#
#     Parameters
#     ----------
#     dfs : dict
#         Dictionary of {key: DataFrame}.
#
#     Returns
#     -------
#     pd.DataFrame
#         Single DataFrame with a MultiIndex where level 0 is the dict key.
#     """
#     return pd.concat(dfs, axis=0)


def join_dict_dataframes(dfs: dict | list, raise_on_conflict: bool = False) -> pd.DataFrame:
    """
    Concatenate DataFrames without adding an extra index level.
    Checks for duplicate index entries across DataFrames before merging.

    Parameters
    ----------
    dfs : dict or list
        Dictionary of {key: DataFrame} or list of DataFrames.
    raise_on_conflict : bool
        If True, raise an error on duplicate indices. If False, warn and keep first occurrence.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame with the original index preserved.
    """
    if not raise_on_conflict:
        return pd.concat(dfs, axis=0)

    else:
        frames = list(dfs.values()) if isinstance(dfs, dict) else dfs

        # Collect all indices and find duplicates
        all_indices = pd.Index([idx for df in frames for idx in df.index])
        duplicates = all_indices[all_indices.duplicated()]

        if not duplicates.empty:
            msg = f"Duplicate index entries found across DataFrames: {duplicates.tolist()}"
            if raise_on_conflict:
                raise ValueError(msg)
            else:
                import warnings
                warnings.warn(msg + " — keeping first occurrence.")

        return pd.concat(frames, axis=0, verify_integrity=raise_on_conflict)


lime.theme.set_style('dark')

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
alpha_folder = project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'
metals_folder = project_folder/'LyC_leakers_COS'/'metals_line_fitting'
lyman_alpha_folder = project_folder/'LyC_leakers_COS'/'Lyman_alpha_fitting'
lsf_obj_folder = project_folder/'LyC_leakers_COS'/'lsf_objects'
results_folder =  project_folder/'LyC_leakers_COS'/'line_measurements'

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

survey_dict = {}
for i, obj in enumerate(target_list):

    if i >= 0:

        # Object folders the files
        print(f'{i}) Galaxy {obj}, z = {cfg_sample['Galaxy_redshifts'][obj]}')
        input_folder_single = obs_folder / 'LyC_leakers_COS' / 'objects_x1d' / f'{obj}'
        output_folder_single = obs_folder / 'LyC_leakers_COS' / 'obj_hasp' / f'{obj}'
        opacityLymanAlpha_df_path = lyman_alpha_folder/f'{obj}_LyAlpha_lines_frame.txt'
        opacity_df_path = metals_folder/f'{obj}_metals_lines_frame.txt'
        voigfitz_reg = metals_folder/f'{obj}_metals_best_fit.reg'
        obj_lsf_file = f"{lsf_obj_folder}/{obj}_hasp_lsf.txt"

        HI_cfg = lyAlpha_data[f'{obj}_LyAlpha_results']
        metals_cfg = cfg_sample['voigtfit_metals_params'][obj]
        alpha_cfg = cfg_sample['voigtfit_lyAlpha_params'][obj]

        # Index the files
        idcs_out = (sample_df.object == obj) & (sample_df.index.get_level_values('state') == 'aspec_manual')

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
        spec.save_spectrum(fname=metals_folder/f'{obj}_metals_spec.txt')

        # Sort the line selection line lists
        list_science, list_masked = unpack_lines(spec, metals_cfg, science_groups=['ISM'],
                                                 mask_groups=['airglow', 'MW', 'nebular', 'stellar'])

        # # Get line list extra lines
        # bands_target = spec.retrieve.lines_frame(band_vsigma=70, **metals_cfg['LyC_lines'])
        #
        # # Plot inputs
        # fig_cfg = {"legend.fontsize": 10,}
        # spec.plot.spectrum(bands=bands_target, line_list=list_masked + list_science, rest_frame=False, in_fig=None, fig_cfg=fig_cfg, ax_cfg={'title': obj})
        # add_opacity_profile(spec, opacity_df_path, voigtfit_pname=voigfitz_reg)

        # Lyman alpha data
        fname = lyman_alpha_folder/f'{obj}_LyAlpha_lines_frame.txt'
        obj_dict = {'alpha': lime.load_frame(fname)}

        for ext in ['metals', 'carbon', 'sulfur']:
            fname = metals_folder/f'{obj}_{ext}_lines_frame.txt'
            if fname.exists():
                print(f'EXTENSION {ext}')
                obj_dict[ext] = lime.load_frame(fname)

        join_df = join_dict_dataframes(obj_dict, raise_on_conflict=True)
        survey_dict[obj] = join_df

LyA_frame = join_dict_dataframes(survey_dict)
LyA_frame.index = LyA_frame.index.set_names(["object", "line"])

LyA_frame.index = LyA_frame.index.set_names(["object", "line"])
LyA_frame = LyA_frame.rename(columns={'ion': 'particle'})

idcs_fine_cIIa = LyA_frame.index.get_level_values('line').isin(['C2_1336A', 'C2_1335.7A'])
LyA_frame.loc[idcs_fine_cIIa, 'particle'] = 'C2*'

version = 'v1'
lime.save_frame(results_folder/f'LyC_voigtprofiles_{version}.txt', LyA_frame)