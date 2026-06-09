import urllib
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt
import toml
import re

from hasp import wrapper
from lmfit.models import GaussianModel, LinearModel
from astropy.coordinates import SkyCoord
from regions import CircleSkyRegion

import shutil
from scipy.constants import m_e
from scipy.signal import fftconvolve
from scipy.signal.windows import gaussian as scipy_gaussian
import lime
from lime import Line
from lime.fitting.lines import c_KMpS
import warnings

vf_database_path = '/home/vital/Dropbox/Astrophysics/Tools/LiMe/LinesDatabase/linelist_voigtfit.dat'
vfdb_df = pd.read_csv(vf_database_path, names=['label', 'ion', 'wl', 'f_trans', 'gam_trans', 'mass_u'], index_col=0,
                                skiprows=5,  sep=r'\s+', skip_blank_lines=True, comment='#')

e_cgs = 4.8032e-10
c_cgs = 2.99792458e10  # cm/s
m_e_grams = m_e * 1000 # grams

DEFAULT_b = {"HI": 10,
             "AlII": 10,
             "CII": 10,
             "CIIa": 10,
             "FeII": 10,
             "NI": 10,
             "NiII": 10,
             "OI": 10,
             "PII": 10,
             "SII": 10,
             "SiII": 10}

DEFAULT_logN = {"HI": 20,
                "AlII": 13,
                "CII": 14.5,
                "CIIa": 14,
                "FeII": 15,
                "NI": 15,
                "NiII": 14,
                "OI": 15,
                "PII": 14,
                "SII": 15,
                "SiII": 14}


def voigt_hjerting_approx(a, x_arr):

    P = x_arr * x_arr
    H_0 = np.exp(-P)
    Q = 1.5/P

    return H_0 - (a / np.sqrt(np.pi) / P) * (H_0 * H_0 * (4 * P * P + 7 * P + 4 + Q) - Q - 1)


def optical_depth_profile(wave_arr, lambda_trans, f_trans, gam_trans, N, b, redshift=0):

    # Convert inputs to CGS if they are in km/s
    b_cgs = b * 1e5
    l0_cm = lambda_trans * 1e-8

    # 1. Cross-section at line center (C_a)
    cs_profile = (np.sqrt(np.pi) * np.square(e_cgs) * f_trans * l0_cm) / (m_e_grams * c_cgs * b_cgs)

    # 2. Damping constant (a)
    damp_const = (l0_cm * gam_trans) / (4 * np.pi * b_cgs)

    # 3. Frequency displacement (x)
    wave_rest = wave_arr / (1 + redshift)
    dop_displ = (c_cgs / b_cgs) * (1.0 - lambda_trans / wave_rest)

    # 4. Final Optical Depth
    tau_arr = cs_profile * N * voigt_hjerting_approx(damp_const, dop_displ)

    return tau_arr



def instrumental_broadening(input_flux, wave_arr, kernel_in, sampling=3, kernel_nsub=6):

    # Adjust the kerneal to the line range
    if isinstance(kernel_in, float):
        dx = np.mean(np.diff(wave_arr))
        xmin = np.log10(wave_arr.min() - 50 * dx)
        xmax = np.log10(wave_arr.max() + 50 * dx)
        N = int(sampling * wave_arr.size)
        profile_wl = np.logspace(xmin, xmax, N)
        pxs = np.diff(profile_wl)[0] / profile_wl[0] * 299792.458
        kernel = kernel_in / pxs / 2.35482

    elif isinstance(kernel_in, np.ndarray):
        N = int(kernel_nsub * len(wave_arr))
        assert kernel_in.shape[0] == N
        # evaluate on the input grid subsampled by `nsub`:
        if kernel_nsub > 1:
            profile_wl = np.linspace(kernel_in.min(), kernel_in.max(), N)
        else:
            profile_wl = kernel_in.copy()

    else:
        err_msg = "Invalid type of `kernel`: %r" % type(kernel_in)
        raise TypeError(err_msg)

    # Kernel from the LSF
    LSF = scipy_gaussian(10*int(kernel) + 1, kernel)
    LSF = LSF / LSF.sum()

    if isinstance(kernel, float):
        profile_broad = fftconvolve(input_flux, LSF, 'same')
        profile_obs = np.interp(x, profile_wl, profile_broad)

    return profile_int


def profile_normflux(tau_arr, wave_arr, kernel, pxs):

    profile_int = np.exp(-tau_arr)

    pad_width = kernel
    padded = np.pad(profile_int, pad_width, mode='edge')  # extend edge values, not zeros

    # Kernel from the LSF
    sigma_instrumental = 1
    LSF = scipy_gaussian(wave_arr.size // 2, sigma_instrumental)
    LSF = LSF / LSF.sum()

    profile_int = fftconvolve(padded, LSF, 'same')
    profile_int = profile_int[pad_width:-pad_width]  #

    # # Perform the convolution
    # fig, ax = plt.subplots()
    # ax.step(wave_arr, profile_int)
    # plt.show()

    return profile_int


def absorption_spectrum(opacity_pname, spec, kernel=20, pxs=1.570539778139733):

    # Check if true
    if not opacity_pname.is_file():
        IOError(f'- WARNING: No opacity lines frame at {opacity_pname}')

    # Container for the opacity
    tau_arr = np.zeros(spec.wave.data.size)

    # Loop though the lines and add up the opacities
    lines_df = lime.load_frame(opacity_pname)
    for row in lines_df.itertuples():
        tau_arr += optical_depth_profile(spec.wave.data,
                                         row.wavelength,
                                         row.osc_str,
                                         row.trans_prob,
                                         np.power(10, row.logN),
                                         row.b,
                                         row.z_line)

    # Compute the normalized flux
    norm_flux = profile_normflux(tau_arr, spec.wave.data, kernel=kernel, pxs=pxs)

    return norm_flux


def lines_frame_to_lime_lines(df, fit_cfg):

    line_list = []
    print(f'\n- Input lines:')
    for label in df.index:
        line = Line.from_transition(label, data_frame=df, fit_cfg=fit_cfg)
        line_list.append(line)
        msg = f' {line.label} {'' if len(line.list_comps) == 1 else line.list_comps}'

    return line_list


def build_ion_dicts(df, param_hdrs=('b', 'z', 'logN', 'var_z', 'var_b', 'var_N', 'tie_z', 'tie_b'), default_dict=None):

    required = {"ion", "origin"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"DataFrame is missing columns: {missing}")

    df = df.copy()

    def _extract_params(group, ion, origin, kinem):
        out = {}
        for col in param_hdrs:
            unique_vals = group[col].dropna().unique()
            if len(unique_vals) == 0:
                out[col] = default_dict[col] if col.startswith('var') or col.startswith('tie') else default_dict[col][ion]
            elif len(unique_vals) == 1:
                out[col] = unique_vals[0]
            else:
                warnings.warn(
                    f"[ion={ion!r}, origin={origin!r}, kinem={kinem}] column {col!r} has "
                    f"{len(unique_vals)} unique values {list(unique_vals)}. "
                    f"Using first: {unique_vals[0]!r}",
                    UserWarning, stacklevel=3,
                )
                out[col] = unique_vals[0]
        return out

    result = {}

    for ion, ion_group in df.groupby("ion", sort=False):
        result[ion] = {}

        for origin, origin_group in ion_group.groupby("origin", sort=False):
            origin_dict = {}

            for kinem, kinem_group in origin_group.groupby("kinem", sort=True):
                params = _extract_params(kinem_group, ion, origin, kinem)
                origin_dict[f"k-{kinem}"] = {'lines_group': kinem_group.index.to_numpy(), 'params': params}

            result[ion][origin] = origin_dict

    return result


# def build_ion_dicts(df, param_hdrs =('b', 'z', 'logN', 'var_z', 'var_b', 'var_N'), default_dict=None):
#     """
#     Return { ion -> { origin -> { col: value } } } for all ions and
#     their observed origins.
#
#     Rules per (ion, origin, column):
#       - One unique non-NaN value  → use it.
#       - Multiple unique non-NaN values → warn, use the first.
#       - All NaN → fall back to DEFAULTS[col].
#
#     Parameters
#     ----------
#     df : pd.DataFrame
#         Must contain columns: "ion", "origin", + all keys in DEFAULTS.
#
#     Returns
#     -------
#     dict  { ion_str: { origin_str: { col: value } } }
#     """
#     required = {"ion", "origin"} | set()
#     missing = required - set(df.columns)
#     if missing:
#         raise ValueError(f"DataFrame is missing columns: {missing}")
#
#     result: dict[str, dict[str, dict]] = {}
#
#     for ion, ion_group in df.groupby("ion", sort=False):
#         result[ion] = {}
#
#         for origin, group in ion_group.groupby("origin", sort=False):
#             origin_dict = {}
#
#             for col in param_hdrs:
#                 unique_vals = group[col].dropna().unique()
#
#                 if len(unique_vals) == 0:
#                     if col.startswith('var'):
#                         origin_dict[col] = default_dict[col]
#                     else:
#                         origin_dict[col] = default_dict[col][ion]
#
#                 elif len(unique_vals) == 1:
#                     origin_dict[col] = unique_vals[0]
#
#                 else:
#                     warnings.warn(f"[ion={ion!r}, origin={origin!r}] column {col!r} has "
#                                             f"{len(unique_vals)} unique values {list(unique_vals)}. "
#                                             f"Using first: {unique_vals[0]!r}",
#                                             UserWarning,
#                                             stacklevel=2,)
#                     origin_dict[col] = unique_vals[0]
#
#             result[ion][origin] = origin_dict
#
#     return result

def _extract_params(section: dict, KEY_PARAMS=['z', 'logN', 'b', 'var_b', 'var_z', 'var_N']) -> dict:
    """Extract only the recognised fit parameters from a config section."""
    return {k: v for k, v in section.items() if k in KEY_PARAMS}



def run_VoigtFit(fpath, spec, lines_df, fit_cfg, conv_dict, resolution=20, obj_redshift=0, output_toml=None,
                 voigt_default_params=None, lsf_file=None, show_plots=False):

    import VoigtFit as vf

    for tie_param in ['tie_b', 'tie_z']:
        if tie_param not in voigt_default_params:
            voigt_default_params[tie_param] = None

    # Generate spectrum object
    obj_name = fpath.stem
    dataset = vf.DataSet(redshift=obj_redshift, name=obj_name)
    dataset.verbose = True

    dataset.add_data(spec.wave.data, spec.flux.data,
                     res=resolution if lsf_file is None else lsf_file,
                     err=spec.err_flux.data,
                     mask=~spec.flux.mask,
                     normalized=True,
                     nsub=6)

    # Loop through the lines, add the regions to the data set and check for the fitting transition configuration
    comps_params = {}
    missing_data = {}
    print('\nInput lines:')
    for line in lines_df.index:
        line = lime.Line.from_transition(line, fit_cfg, lines_df)
        print(f'\n{line}' + f' {line.list_comps if len(line.list_comps) > 1 else ""}')

        # Compute the kernel for the function
        for trans in line.list_comps:

            # Compile the lines missing the atomic data
            if (trans.atom_data.osc_str is None) or (trans.atom_data.trans_prob) is None:
                missing_data[label_vf] = trans.label

            # Identifiers
            label_vf = conv_dict[trans.core]
            ion = trans.classic_notation(just_particle=True)
            origin = 'none' if trans.origin is None else trans.origin
            kinem = trans.kinem
            z = spec.redshift if trans.redshift in [None, 'none'] else trans.redshift

            # Local values
            b = fit_cfg.get('b', {}).get(ion)
            logN = fit_cfg.get('logN', {}).get(ion)

            # Add line with its region
            if label_vf not in dataset.all_lines:
                v_blue = c_KMpS * (line.wavelength - line.mask[2]) / line.wavelength
                v_red = c_KMpS * (line.mask[3] - line.wavelength) / line.wavelength
                dataset.add_line(line_tag=label_vf, velspan=(-v_blue, v_red))
                print(f'- {trans} = {label_vf} ({-v_blue:0.1f} km/s, {v_red:0.1f} km/s)')

            # Get tied kinematics
            parent = fit_cfg.get('kinem', {}).get(origin, {}).get(f'k-{kinem}_tie')
            if parent is not None:
                parent_line = lime.Line.from_transition(parent)
                if parent_line.classic_notation(just_particle=True) == ion:
                    parent = None
            tie_z, tie_b = parent, parent

            # Container with the line data
            local_params = dict(vf_label=label_vf,
                                vf_idx = None,
                                ion=ion,
                                wavelength = trans.wavelength,
                                origin = origin,
                                kinem = kinem,
                                z = z,
                                b = voigt_default_params['b'][ion] if b is None else b,
                                logN = voigt_default_params['logN'][ion] if logN is None else logN,
                                var_z = True,
                                var_b = True,
                                var_N = True,
                                tie_z = tie_z,
                                tie_b = tie_b,
                                idx_comp = np.nan,
                                group_label = trans.group_label if trans.group_label else 'none',
                                line = line,
                                ion_lime = trans.particle.label,
                                osc_str = trans.atom_data.osc_str,
                                trans_prob = trans.atom_data.trans_prob)

            # Global configuration (Priority to update Origin < Kinematic < Ion)
            grouped_params = {}
            if 'kinem' in fit_cfg:

                if origin in fit_cfg['kinem']:

                    # Origin params
                    grouped_params.update(_extract_params(fit_cfg['kinem'][origin]))

                    # Kinematic group
                    kinem_key = f'k-{kinem}'
                    if kinem_key in fit_cfg['kinem'][origin]:
                        grouped_params.update(fit_cfg['kinem'][origin][kinem_key])

                    if ion in fit_cfg['kinem'][origin]:
                        grouped_params.update(fit_cfg['kinem'][origin][ion])


            # Grouped parameters update transition data (in voigtfit)
            local_params.update(grouped_params)

            # Store the parameters and
            comps_params[trans] = local_params.copy()

            # Create an identifier
            identifier = f'{comps_params[trans]["ion"]}_{comps_params[trans]["kinem"]}_{comps_params[trans]["origin"]}_{comps_params[trans]["z"]}'
            comps_params[trans]['lime_label'] = identifier

    # Stop for missing data
    if len(missing_data) > 0:
        print(f'\n- MISSING ATOMIC DATA:')
        for key, value in missing_data.items():
            print(f'{key} -> {value}')
        raise KeyError(f'Missing atomic data')

    # Warn if missig atomic data
    if (trans.atom_data.osc_str is None) or (trans.atom_data.trans_prob) is None:
        missing_data[label_vf] = trans.label

    # Determine components and resolve conflicts
    df_lines = pd.DataFrame(index=comps_params.keys(), data=comps_params.values())
    # df_lines['vf_idx'] = df_lines.groupby(['origin', 'kinem']).ngroup()
    df_lines['vf_idx'] = pd.factorize(pd.MultiIndex.from_arrays([df_lines['origin'], df_lines['kinem']]))[0]

    # Cleanup the tie_z
    idcs_tie = pd.notnull(df_lines.tie_z)
    for idx in df_lines.loc[idcs_tie].index:
        parent = df_lines.loc[idx, 'tie_z']
        # parent_ion, parent_order = df_lines.loc[parent, ['ion', 'kinem']]
        parent_ion, parent_order = df_lines.loc[parent, ['ion', 'vf_idx']]
        df_lines.loc[idx, 'tie_z'] = f'z{parent_order}_{parent_ion}'
        df_lines.loc[idx, 'tie_b'] = f'b{parent_order}_{parent_ion}'

    comps_dict = build_ion_dicts(df_lines, default_dict=voigt_default_params)
    print(df_lines.to_string())

    # Check for conflicts and components
    print(f'\n- Components:')
    for ion in comps_dict.keys():
        # print(f'-- {ion}: ')
        idx_comp = 0
        for origin, kinem_comps in comps_dict[ion].items():
            for kinem_label, comp_cfg in kinem_comps.items():
                lines_group = comp_cfg['lines_group']
                params = comp_cfg['params']
                # print(f'--- Component {idx_comp}) {lines_group}: ',", ".join(f"{k}: {v}" for k, v in params.items()))
                df_lines.loc[df_lines.index.isin(lines_group), 'idx_comp'] = idx_comp
                df_lines.loc[df_lines.index.isin(lines_group), 'z'] = params['z']
                dataset.add_component(ion=ion, **params)
                idx_comp += 1

    # Missing message:
    if len(missing_data) > 0:
        print(f'\n- MISSING ATOMIC DATA:')
        for key, value in missing_data.items():
            print(f'{key} -> {value}')

    for ion, comp in dataset.components.items():
        print(ion)
        for kine_c in comp:
            print(f'-- {kine_c}, var_z={kine_c.var_z}, var_b={kine_c.var_b}, var_N={kine_c.var_N}, '
                  f'             tie_z={kine_c.tie_z}, tie_b={kine_c.tie_b}, tie_N={kine_c.tie_N},')

    # Normalize and mask if necessary
    dataset.prepare_dataset(norm=False, mask=False)

    # Fit the data
    _, _ = dataset.fit(verbose=True)

    # Save measurements in LiMe results
    frame_path = fpath.parent/f'{obj_name}_lines_frame.txt'
    extract_voigtfit_dataframe(frame_path, df_lines, comps_dict, dataset)

    # LiMe save measurements
    param_dict = extract_voigtfit_params(obj_name, df_lines, comps_dict, dataset)
    lime.save_cfg(output_toml, param_dict, section_name=f'{obj_name}_results', clear_section=True)
    strip_all_quotes(output_toml)

    # Voigtfit measurements
    dataset.save_fit_regions(filename=str(fpath.parent/f'{obj_name}_best_fit'))

    # Save the results
    if show_plots:
        dataset.plot_fit(individual=True, xunit='wl')
        plt.show(block=True)
    else:
        dataset.plot_fit(filename=str(fpath.parent / f'{obj_name}_voigtfit_plot.pdf'), individual=False, xunit='vel')

    # -- Print total column densities
    dataset.print_total()
    dataset.save_parameters(str(fpath.parent/f'{obj_name}_voigtfit_parameters'))
    # dataset.save_cont_parameters_to_file(str(fpath.parent/f'{obj_name}_voigtfit_cont'))

    return


# def run_VoigtFit(fpath, spec, lines_df, fit_cfg, conv_dict, resolution=20, obj_redshift=0, output_toml=None,
#                  voigt_default_params=None, lsf_file=None, show_plots=False):
#
#     import VoigtFit as vf
#
#     # Generate spectrum object
#     obj_name = fpath.stem
#     dataset = vf.DataSet(redshift=obj_redshift, name=obj_name)
#     dataset.verbose = True
#
#     dataset.add_data(spec.wave.data, spec.flux.data,
#                      res=resolution if lsf_file is None else lsf_file,
#                      err=spec.err_flux.data,
#                      mask=~spec.flux.mask,
#                      normalized=True,
#                      nsub=6)
#
#     # Loop through the lines, add the regions to the data set and check for the fitting transition configuration
#     comps_params = {}
#     missing_data = {}
#     print('\nInput lines:')
#     for line in lines_df.index:
#         line = lime.Line.from_transition(line, fit_cfg, lines_df)
#         msg = f'\n{line}' + f' {line.list_comps if len(line.list_comps) > 1 else ""}'
#         print(msg)
#
#         # Compute the kernel for the function
#         for trans in line.list_comps:
#             label_vf = conv_dict[trans.core]
#
#             # Add line with its region
#             if label_vf not in dataset.all_lines:
#                 v_blue = c_KMpS * (line.wavelength - line.mask[2]) / line.wavelength
#                 v_red = c_KMpS * (line.mask[3] - line.wavelength) / line.wavelength
#                 print(f'- {trans} = {label_vf} ({-v_blue:0.1f} km/s, {v_red:0.1f} km/s)')
#                 dataset.add_line(line_tag=label_vf, velspan=(-v_blue, v_red))
#
#             # Kinematic parameters
#             kinem = trans.kinem
#
#             # Origin parameters
#             orig = 'none' if trans.origin is None else trans.origin
#             z_comp = fit_cfg.get(trans, {}).get('z')
#             if z_comp is None:
#                 z_comp = spec.redshift if trans.redshift in [None, 'none'] else trans.redshift
#
#             # Ion based paramters
#             ion = trans.classic_notation(just_particle=True)
#
#             # Transition based parameters...
#             b = fit_cfg.get(trans, {}).get('b', np.nan)
#             logN = fit_cfg.get(trans, {}).get('logN', np.nan)
#             var_z_i = fit_cfg.get(trans, {}).get('var_z')
#             var_b_i = fit_cfg.get(trans, {}).get('var_b')
#             var_N_i = fit_cfg.get(trans, {}).get('var_N')
#
#             # Ion based parameters
#             if np.isnan(b) and 'b' in fit_cfg.get(trans.particle.label, {}):
#                 b = fit_cfg[trans.particle.label]['b']
#             if np.isnan(logN) and 'logN' in fit_cfg.get(trans.particle.label, {}):
#                 logN = fit_cfg[trans.particle.label]['logN']
#
#             comp_label = f'{ion}_{kinem}_{orig}_{z_comp}'
#             comps_params[trans] = dict(vf_label=label_vf,
#                                        idx_comp=np.nan,
#                                        wavelength=trans.wavelength,
#                                        lime_label=comp_label,
#                                        group_label=trans.group_label if trans.group_label else 'none',
#                                        line=line,
#                                        origin=orig,
#                                        kinem=kinem,
#                                        ion=ion,
#                                        ion_lime=trans.particle.label,
#                                        b=b,
#                                        z=np.nan if z_comp is None else z_comp,
#                                        logN=logN,
#                                        var_z=var_z_i,
#                                        var_b=var_b_i,
#                                        var_N=var_N_i,
#                                        osc_str=trans.atom_data.osc_str,
#                                        trans_prob=trans.atom_data.trans_prob,
#                                        )
#
#             # Warn if missig atomic data
#             if (trans.atom_data.osc_str is None) or (trans.atom_data.trans_prob) is None:
#                 missing_data[label_vf] = trans.label
#
#     # Prepare default values
#     voigt_default_params = {**voigt_default_params}
#     voigt_default_params['var_z'] = fit_cfg.get('var_z', True)
#     voigt_default_params['var_b'] = fit_cfg.get('var_b', True)
#     voigt_default_params['var_N'] = fit_cfg.get('var_N', True)
#
#     # Determine components and resolve conflicts
#     df_lines = pd.DataFrame(index=comps_params.keys(), data=comps_params.values())
#     comps_dict = build_ion_dicts(df_lines, default_dict=voigt_default_params)
#     for ion, ion_conf in comps_dict.items():
#         print(f'Ion: {ion}')
#         for origin, origin_conf in ion_conf.items():
#             print(f'- Origin: {origin}')
#             for kinem, k in origin_conf.items():
#                 print(f'-- {kinem}: {k}')
#
#     # Check for conflicts and components
#     print(f'- Components:')
#     for ion in comps_dict.keys():
#         print(f'-- {ion}: ')
#         idx_comp = 0
#         for origin, kinem_comps in comps_dict[ion].items():
#             for kinem_label, comp_cfg in kinem_comps.items():
#                 lines_group = comp_cfg['lines_group']
#                 params = comp_cfg['params']
#                 print(f'--- Component {idx_comp}) {lines_group}: ',", ".join(f"{k}: {v}" for k, v in params.items()))
#                 df_lines.loc[df_lines.index.isin(lines_group), 'idx_comp'] = idx_comp
#                 dataset.add_component(ion=ion, **params)
#                 idx_comp += 1
#
#     # Add tied up components
#     if 'kinem' in fit_cfg:
#         if 'ion' in fit_cfg['kinem']:  # Tied by species:
#             for ion, tied_species in fit_cfg['kinem']['ion'].items():
#                 for child_ion in tied_species:
#                     if dataset.has_ion(ion):
#                         print(f'Replacing components for {child_ion} with those from {ion}')
#                     dataset.copy_components(child_ion, ion, logN=voigt_default_params['logN'][child_ion], tie_z=True, tie_b=False)
#
#
#
#     print('\nInput Voigtfit lines')
#     print(dataset.lines)
#
#     print('\nInput Voigtfit components')
#     for ion, comp_list in dataset.components.items():
#         print(f'Ion {ion}')
#         for comp in comp_list:
#             print(f'-- z={comp.z} var_z={comp.var_z}; z={comp.b} var_b={comp.var_b}; logN={comp.logN} var_N={comp.var_N}')
#
#     # Missing message:
#     if len(missing_data) > 0:
#         print(f'\n- MISSING ATOMIC DATA:')
#         for key, value in missing_data.items():
#             print(f'{key} -> {value}')
#
#     # Normalize and mask if necessary
#     dataset.prepare_dataset(norm=False, mask=False)
#
#     # Fit the data
#     popt, chi2 = dataset.fit(verbose=True, factor=10.)
#
#     # Save measurements in LiMe results
#     frame_path = fpath.parent/f'{obj_name}_lines_frame.txt'
#     extract_voigtfit_dataframe(frame_path, df_lines, comps_dict, dataset)
#
#     # LiMe save measurements
#     param_dict = extract_voigtfit_params(obj_name, df_lines, comps_dict, dataset)
#     lime.save_cfg(output_toml, param_dict, section_name=f'{obj_name}_results', clear_section=True)
#     strip_all_quotes(output_toml)
#
#     # Voigtfit measurements
#     dataset.save_fit_regions(filename=str(fpath.parent/f'{obj_name}_best_fit'))
#
#     # Save the results
#     if show_plots:
#         dataset.plot_fit(individual=True, xunit='wl')
#         plt.show(block=True)
#     else:
#         dataset.plot_fit(filename=str(fpath.parent / f'{obj_name}_voigtfit_plot.pdf'), individual=False, xunit='wl')
#
#     # -- Print total column densities
#     dataset.print_total()
#     dataset.save_parameters(str(fpath.parent/f'{obj_name}_voigtfit_parameters'))
#     # dataset.save_cont_parameters_to_file(str(fpath.parent/f'{obj_name}_voigtfit_cont'))
#
#     return


def plot_profile_components(output_fpath, object_label, spec, dataset):


    cfg_fig = {"figure.figsize": (10, 5), "figure.dpi": 350}

    spec.plot.spectrum(in_fig=None, show_cont=True, fig_cfg=cfg_fig)
    spec.plot.ax.set_xlim(spec.wave.min(), spec.wave.max())
    spec.plot.ax.set_ylim(-0.1, 1.4)

    tau_arr = np.zeros(spec.wave.data.size)
    for j, ion in enumerate(sorted(dataset.components.keys())):
        for i in range(len(dataset.components[ion])):
            z = dataset.best_fit[f'z{i}_{ion}']
            logN = dataset.best_fit[f'logN{i}_{ion}']
            b = dataset.best_fit[f'b{i}_{ion}']
            f_trans, gam_trans = vfdb_df.loc['HI_1215', ['f_trans', 'gam_trans']]

            tau_arr += optical_depth_approx(spec.wave.data, 1215.6696, f_trans, gam_trans, 10 ** logN, b.value, z.value)
            if i == 0:
                try:
                    label = f'{object_label.split('_')[0]}, z={z.value:0.3f}±{logN.stderr:0.3f}, logN={logN.value:0.3f}±{logN.stderr:0.3f}, b={b.value:0.3f}±{b.stderr:0.3f}'
                except:
                    label = f'{object_label.split('_')[0]}'
            else:
                try:
                    label = f'Milky way, z={z.value:0.3f}±{z.stderr:0.3f}, logN={logN.value:0.3f}±{logN.stderr:0.3f}, b={b.value:0.3f}±{b.stderr:0.3f}'
                except:
                    label = 'Milky way'

    norm_profile = instrumental_broadening(tau_arr, spec.wave.data, kernel=20, pxs=1.570539778139733)
    spec.plot.ax.plot(spec.wave.data, norm_profile, linestyle='--', linewidth=2, label=label)

    spec.plot.ax.legend(loc=9)
    plt.tight_layout()
    plt.savefig(output_fpath)
    np.savetxt(Path(output_fpath).with_suffix('.txt'), np.column_stack([spec.wave.data, norm_profile]))


    return


def extract_voigtfit_params(obj_label, inputs_df, input_comps, dataset, velocity=False):

    """
    Convert best-fit parameters to a single dictionary section
    named after the variable prefix.
    """

    z_sys = dataset.redshift

    # Everything goes into this one dictionary
    combined_data = {"z_sys": float(z_sys), "chi2": float(dataset.chi2), "n_free": int(dataset.minimizer.result.nfree)}

    # Add the components to the same dictionary
    for ion in sorted(dataset.components.keys()):
        for i in range(len(dataset.components[ion])):
            best_fit = dataset.best_fit
            z = dataset.best_fit[f'z{i}_{ion}']
            logN = dataset.best_fit[f'logN{i}_{ion}']
            b = dataset.best_fit[f'b{i}_{ion}']
            lambda_i = dataset.best_fit[f'b{i}_{ion}']

            # Handle Redshift vs Velocity
            if velocity:
                val = (z.value - z_sys) / (z_sys + 1) * 299792.458
                err = (z.stderr / (z_sys + 1) * 299792.458) if z.stderr else 0.0
                param_name = "velocity_kms"
            else:
                val = z.value
                err = z.stderr if z.stderr else 0.0
                param_name = "redshift"

            # # Dot notation keys added directly to the main dict
            # combined_data[f"components.{ion}.{i}.{param_name}"] = [float(val), float(err)]
            # combined_data[f"components.{ion}.{i}.b_kms"] = [float(b.value), float(b.stderr) if b.stderr else 0.0]
            # combined_data[f"components.{ion}.{i}.logN"] = [float(logN.value), float(logN.stderr) if logN.stderr else 0.0]


    # Dataframe with the results

    # Return as a nested dict so TOML creates one header: [prefix]
    return combined_data


def extract_voigtfit_dataframe(fname, inputs_df, input_comps, dataset, velocity=False):

    """
    Convert best-fit parameters to a single dictionary section
    named after the variable prefix.
    """
    # Results container
    results_dict = {}

    # Global parameters
    chi2 = dataset.chi2
    z_sys = dataset.redshift
    n_free = dataset.minimizer.result.nfree
    success_check = dataset.minimizer.result.success

    # Add the components to the same dictionary
    for ion in sorted(dataset.components.keys()):
        results_dict[ion] = {}
        for i in range(len(dataset.components[ion])):
            z = dataset.best_fit[f'z{i}_{ion}']
            logN = dataset.best_fit[f'logN{i}_{ion}']
            b = dataset.best_fit[f'b{i}_{ion}']

            # Handle Redshift vs Velocity
            v_r = (z.value - z_sys) / (z_sys + 1) * 299792.458
            v_r_err = (z.stderr / (z_sys + 1) * 299792.458) if z.stderr else np.nan

            results_dict[ion][i] = {}
            results_dict[ion][i]['b'] = b.value
            results_dict[ion][i]['b_err'] = float(b.stderr) if b.stderr else np.nan

            results_dict[ion][i]['logN'] = logN.value
            results_dict[ion][i]['logN_err'] = float(logN.stderr) if logN.stderr else np.nan

            results_dict[ion][i]['z_line'] = z.value
            results_dict[ion][i]['z_line_err'] = z.stderr if z.stderr else np.nan

            results_dict[ion][i]['v_r'] = v_r
            results_dict[ion][i]['v_r_err'] = v_r_err

    # Parse the original lines to their components
    columns = ['wavelength', 'ion', 'idx_comp', 'osc_str', 'trans_prob', 'group_label',
               'b', 'b_err', 'logN', 'logN_err', 'v_r', 'v_r_err', 'z_line', 'z_line_err']
    common_hdrs = ['wavelength', 'osc_str', 'trans_prob', 'group_label']
    results_hdrs = ['b', 'b_err', 'logN', 'logN_err', 'v_r', 'v_r_err', 'z_line', 'z_line_err']

    out_df = pd.DataFrame(columns=columns)
    for line in inputs_df.index:
        ion, origin = inputs_df.loc[line, ['ion', 'origin']]
        idx_comp = int(inputs_df.loc[line, 'idx_comp'])
        idcs_rows = (inputs_df.ion == ion) & (inputs_df.origin == origin)

        out_df.loc[line, common_hdrs] = inputs_df.loc[line, common_hdrs]
        out_df.loc[line, 'ion'] = inputs_df.loc[line, 'ion_lime']
        out_df.loc[line, 'idx_comp'] = inputs_df.loc[line, 'idx_comp']

        for hdr in results_hdrs:
            out_df.loc[idcs_rows, hdr] = results_dict[ion][idx_comp][hdr]

    # Global parameters
    out_df['observations'] = 'none' if success_check else 'failed_convergence'
    out_df['chi2'] = chi2

    # Save to a frame file
    lime.save_frame(fname, out_df, **{'z_sys': z_sys})
    print(out_df.to_string())

    return

# def extract_voigtfit_params(dataset, prefix, velocity=False):
#     """
#     Convert best-fit parameters to a dictionary where section names
#     start with a custom variable string prefix.
#     """
#     z_sys = dataset.redshift
#
#     # Define your section names using the variable prefix
#     meta_section = f"{prefix}_results.metadata"
#     comp_section = f"{prefix}_results.components"
#
#     # Root dictionary using dynamic keys
#     results = {meta_section: {"z_sys": float(z_sys),
#                               "chi2": float(dataset.chi2),
#                               "n_free": int(dataset.minimizer.result.nfree)
#                              },
#                comp_section: {}}
#
#     for ion in sorted(dataset.components.keys()):
#         for i in range(len(dataset.components[ion])):
#             z = dataset.best_fit[f'z{i}_{ion}']
#             logN = dataset.best_fit[f'logN{i}_{ion}']
#             b = dataset.best_fit[f'b{i}_{ion}']
#
#             # Handle Redshift vs Velocity
#             if velocity:
#                 val = (z.value - z_sys) / (z_sys + 1) * 299792.458
#                 err = (z.stderr / (z_sys + 1) * 299792.458) if z.stderr else 0.0
#                 param_name = "velocity_kms"
#             else:
#                 val = z.value
#                 err = z.stderr if z.stderr else 0.0
#                 param_name = "redshift"
#
#             # Populate the components with the dot notation inside the prefixed section
#             results[comp_section][f"{ion}.{i}.{param_name}"] = [float(val), float(err)]
#             results[comp_section][f"{ion}.{i}.b_kms"] = [float(b.value), float(b.stderr) if b.stderr else 0.0]
#             results[comp_section][f"{ion}.{i}.logN"] = [float(logN.value), float(logN.stderr) if logN.stderr else 0.0]
#
#     return results
def strip_all_quotes(filename):
    # Read the whole file into memory
    with open(filename, 'r') as f:
        content = f.read()

    # Remove every single double quote character
    clean_content = content.replace('"', '')

    # Save it back to the same file
    with open(filename, 'w') as f:
        f.write(clean_content)

    print(f" - All quotes removed from {filename}")


def clean_toml_section_quotes(filename, section_name):
    """
    Removes quotes from keys only within a specific [section_name].
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    new_lines = []
    inside_target_section = False

    # Regex for a TOML section header: e.g., [my_section]
    section_pattern = re.compile(r'^\[(.+)\]')
    # Regex for quoted keys: "HI.0.redshift" =
    quote_pattern = re.compile(r'^"([\w\.]+)"\s*=')

    for line in lines:
        # Check if we are entering a new section
        section_match = section_pattern.match(line.strip())
        if section_match:
            current_section = section_match.group(1)
            # Toggle 'inside' flag based on the section name
            inside_target_section = (current_section == section_name)

        # If we are in the target section, strip the quotes from keys
        if inside_target_section:
            line = quote_pattern.sub(r'\1 =', line)

        new_lines.append(line)

    with open(filename, 'w') as f:
        f.writelines(new_lines)


def run_VoigtFit_por_si(obj_name, objSpec, line_list, comps_dict, conv_dict, resolution=20):


    import VoigtFit as vf

    # Generate spectrum object
    dataset = vf.DataSet(0, name=obj_name)
    dataset.verbose = True

    dataset.add_data(objSpec.wave.data, objSpec.flux.data,
                     res=resolution,
                     err=objSpec.err_flux.data,
                     mask=~objSpec.flux.mask,
                     normalized=True,
                     nsub=6)

    # Declare the lines
    for i, line in enumerate(line_list):
        v_blue = c_KMpS * (line.wavelength - line.mask[2]) / line.wavelength
        v_red = c_KMpS * (line.mask[3] - line.wavelength) / line.wavelength
        dataset.add_line(line_tag=conv_dict[line.core], velspan=(-v_blue, v_red))
        msg = f' {line.label} {'' if len(line.list_comps) == 1 else line.list_comps} ({-v_blue}, {v_red}) km/s'
        print(msg)

    # Define kinematic components
    if 'kinem_comps' in comps_dict:
        for line_label in comps_dict['kinem_comps'].keys():
            if line_label in line_list:
                line_ref = line_list[line_list.index(line_label)]
                for j, trans in enumerate(line_ref.list_comps):
                    kin_params = {}
                    kin_params['ion'] = trans.particle.classic_notation
                    kin_params['z'] = objSpec.redshift if trans.redshift is None else trans.redshift
                    for i, (param, value_list) in enumerate(comps_dict['kinem_comps'][line_ref.label].items()):
                        kin_params[param] = value_list[j]
                    dataset.add_component(**kin_params)

    # Map kinematic components
    if 'kinem_export' in comps_dict:
        for line_label in comps_dict['kinem_export'].keys():
            if line_label in line_list:
                line_ref = line_list[line_list.index(line_label)]
                export_dict = comps_dict['kinem_export'][line_ref.label]
                parent_ion = line_ref.particle.classic_notation
                for j, child_ion in enumerate(export_dict['to_ion']):
                    dataset.copy_components(from_ion=parent_ion, to_ion=child_ion,
                                            tie_b=export_dict['tie_b'][j], tie_z=export_dict['tie_z'][j])

    # Normalize and mask if necessary
    dataset.prepare_dataset(norm=False, mask=False)

    # Fit the data
    popt, chi2 = dataset.fit(verbose=True, factor=10.)

    # Save the results
    dataset.plot_fit(filename=f'{dataset.name}_voigtfit')

    # -- Print total column densities
    dataset.print_total()

    return


def unpack_lines(spec, obj_cfg, science_groups, mask_groups, plot_lines=False):

    # Loop through the origins and add the lines
    out_dict = {0: None, 1: None}

    for i, current_group_names in enumerate([science_groups, mask_groups]):
        group_lines = []
        for orig in current_group_names:
            if orig in obj_cfg:
                line_list = obj_cfg[orig]['line_list']
                line_list = line_list if 'ignore_list' not  in obj_cfg[orig] else list(set(line_list) - set(obj_cfg[orig]['ignore_list']))
                group_lines += lime.Line.from_list(line_list, obj_cfg[orig].get('fit_cfg'), origin=orig)
        out_dict[i] = group_lines

    bands_science = spec.retrieve.lines_frame(line_list=out_dict[0], **obj_cfg['science_bands'])
    bands_mask = spec.retrieve.lines_frame(line_list=out_dict[1], **obj_cfg['mask_bands'])

    if plot_lines:
        spec.plot.spectrum(line_list=out_dict[0]+out_dict[1])

    list_science = lime.Line.from_list(bands_science.index.to_numpy(), data_frame=bands_science)
    list_masked = lime.Line.from_list(bands_mask.index.to_numpy(), data_frame=bands_mask)

    return list_science, list_masked


### trying scipy instead ####
def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-(((x - mean) / stddev) ** 2) / 2)


### trying the function implemented in IDL gaussfit ####
def fit_func(x, a0, a1, a2, a3, a4, a5):
    z = (x - a1) / a2
    y = a0 * np.exp(-z ** 2 / 2) + a3 + a4 * x + a5 * x ** 2
    return y


def run_hasp_wrapper(in_folder, out_folder, cross_program=True):

    # Run HASP wrapper
    if out_folder.is_dir():
        print(f'\n2) Cross program ', out_folder)
        wrapper.main(in_folder, outdir=out_folder, clobber=True, cross_program=cross_program)
    else:
        raise BrokenPipeError(f'Output dir not found {out_folder}')

    return


def list_files_with_extension(base_dir, extension):
    """
    Return a list of files with a given extension under base_dir, recursively.

    Args:
        base_dir (str or Path): The directory to search.
        extension (str): File extension (e.g. ".txt", ".fits").

    Returns:
        list[Path]: List of Path objects pointing to the files.
    """
    base_path = Path(base_dir)
    return [f for f in base_path.rglob(f"*{extension}") if f.is_file()]


def add_cos_obs(file_list, df, sample_name, state, ref_dict):

    for spec_file in file_list:
        print(spec_file)
        hdr0 = fits.getheader(spec_file, ext=0)
        hdr1 = fits.getheader(spec_file, ext=1)

        row_data = {}
        row_data['sample'] = f'{sample_name}_{hdr0['PROPOSID']}'
        row_data['id'] = f"{hdr0['TARGNAME']}"

        # Extra ID:
        if 'ASN_ID' in hdr0:
            offset_id = f"{hdr0['ASN_ID']}"
        elif 'IPPPSSOO' in hdr0:
            offset_id = f"{hdr0['IPPPSSOO']}_{hdr0['GRATING']}"
        else:
            offset_id = f"acqim_{hdr0['ROOTNAME']}"

        if state == 'x1d':
            offset_id = f"{hdr0['ROOTNAME']}"

        if state == 'counts':
            offset_id = f"{spec_file.parts[-2]}_{hdr0['ROOTNAME']}"

        if state in ['flt', 'cal']:
            offset_id = f"{spec_file.parts[-2]}_{hdr0['ROOTNAME']}"

        if 'MULTI' in row_data['id']:
            offset_id = f"NUM_EXP{hdr0['NUM_EXP']}_PINAME{hdr0['PINAME']}"

        row_data['offset_id'] = offset_id

        row_data['state'] = state

        # row_data['object'] = hdr0['TARGNAME'].lower()
        if hdr0['TARGNAME'] in ref_dict['Galaxy_alias']:
            row_data['object'] = hdr0['TARGNAME']
        else:
            row_data['object'] = next((k for k, v in ref_dict['Galaxy_alias'].items() if hdr0['TARGNAME'] in v), None)

        if row_data['object'] is None:
            print(hdr0['TARGNAME'])

        if 'MULTI' in row_data['id']:
            row_data['object'] = row_data['id'].split('_')[0]

        if row_data['object'] is None and state == 'counts':
            print('MISSINGGGG ', row_data['id'])

        row_data['RA'] = f"{hdr0['RA_TARG']}" if 'RA_TARG' in hdr0 else None
        row_data['DEC'] = f"{hdr0['DEC_TARG']}" if 'RA_TARG' in hdr0 else None
        row_data['PID'] = hdr0['PROPOSID']
        row_data['redshift'] = ref_dict['Galaxy_redshifts'].get(row_data['object'])

        row_data['instr'] =  hdr0['instrume']

        # Grating
        if 'OPT_ELEM' in hdr0:
            row_data['grating'] = f"{hdr0['OPT_ELEM']}"
        elif 'GRATING' in hdr0:
            row_data['grating'] = hdr0['GRATING']
        else:
            row_data['grating'] = f"{hdr0['INSTRUME']}"

        row_data['filecode'] =  hdr0['ASN_ID'] if 'ASN_ID' in hdr0  else spec_file.parts[-2]
        row_data['filepath'] =  Path(*spec_file.parts[5:])

        row_data['cenwave'] = hdr0.get("CENWAVE")
        row_data['life_adj'] = hdr0.get("LIFE_ADJ")
        row_data['disptab'] = hdr0.get("DISPTAB")
        row_data['detector'] = hdr0.get("DETECTOR")

        # Append
        if row_data['grating'] != 'G140L':
            if not isinstance(df.index, pd.MultiIndex):   # True
                df.loc[len(df)] = row_data
            else:
                idx = (row_data["sample"], row_data["id"], row_data["offset_id"], row_data["state"])
                if idx not in df.index:
                    df.loc[idx] = {k: v for k, v in row_data.items() if k not in df.index.names}

    return


def fetch_files(det, grating, lpPos, cenwave, disptab, datadir):


    """
    Given all the inputs, this will download both
    the LSF and Disptab files to use in the convolution and return their paths.

    Input:
    det (str): The detector used
    grating (str): Type of grating used
    lpPos (str): Lifetime position used
    cenwave (str): Central wavelength used
    disptab (str): DISPTAB used (will get the path in the function)

    Returns:
    LSF_file_name (str): filename of the new downloaded LSF file
    disptab_path (str): path to the new downloaded disptab file
    """

    # Link to where all the files live
    COS_site_rootname = ("https://www.stsci.edu/files/live/sites/www/files/"
                        "home/hst/instrumentation/cos/"
                        "performance/spectral-resolution/_documents/")

    # Only one file for NUV
    if det == "NUV":
        LSF_file_name = "nuv_model_lsf.dat"

    # FUV files follow a naming pattern
    elif det == "FUV":
        LSF_file_name = f"aa_LSFTable_{grating}_{cenwave}_LP{lpPos}_cn.dat"

    # Where to find file online
    LSF_file_webpath = COS_site_rootname + LSF_file_name
    urllib.request.urlretrieve(LSF_file_webpath, str(datadir / LSF_file_name))

    # Where to save file to locally
    print(f"Downloaded LSF file to {str(datadir / LSF_file_name)}")

    # # And we'll need to get the DISPTAB file as well
    # disptab_path = str(datadir / disptab)
    # urllib.request.urlretrieve(f"https://hst-crds.stsci.edu/unchecked_get/references/hst/{disptab}", disptab_path)
    #
    # print(f"Downloaded DISPTAB file to {disptab_path}")

    return LSF_file_name


def image_gauss_fitting(y):

    # x grid (replace with your wavelength array if you have one)
    x = np.arange(len(y), dtype=float)

    # ------------------------------
    # Build model: two Gaussians + linear continuum (m*x + b)
    # ------------------------------
    g1 = GaussianModel(prefix='g1_')
    cont = LinearModel(prefix='c_')  # parameters: c_slope, c_intercept
    model = g1 + cont
    params = model.make_params()

    # Rough linear baseline from endpoints (robust to peaks)
    m0 = (np.median(y[-2:]) - np.median(y[:2])) / (x[-1] - x[0] + 1e-12)
    b0 = np.median(y[:2]) - m0 * x[0]
    params['c_slope'].set(m0, vary=False)
    params['c_intercept'].set(b0, vary=False)

    # Centers (bounded to data range)
    idx_max = np.argmax(y)
    span = (x.max() - x.min()) / 8.0
    params['g1_amplitude'].set(value=y[idx_max], min=0)
    params['g1_center'].set(value=idx_max, min=x.min(), max=x.max())
    params['g1_sigma'].set(value=span, min=1)

    # Fit
    result = model.fit(y, params, x=x)

    # Plot
    xx = np.linspace(x.min(), x.max(), 2000)
    y_fit = result.eval(x=xx)
    comp = model.eval_components(params=result.params, x=xx)

    # Results
    arr_dict = {'c_': comp['c_'], 'g1_': comp['g1_'], 'g1_fwhm': result.params['g1_fwhm'].value}

    return arr_dict


def generate_apperture_mask(imshape, wcs, ra, dec, radius, mask_points):

    # Generate the circular mask
    sky_region = CircleSkyRegion(center=SkyCoord(ra=ra, dec=dec, frame="icrs"), radius=radius)
    circle_pix = sky_region.to_pixel(wcs=wcs)
    circle_mask = circle_pix.to_mask(mode="center").to_image(imshape).astype(bool)

    if mask_points is not None:

        # Convert points to coordinates
        c1 = SkyCoord(mask_points[0][0], mask_points[0][1])
        c2 = SkyCoord(mask_points[1][0], mask_points[1][1])

        # Generate the line mask
        yy, xx = np.mgrid[:imshape[0], :imshape[1]]
        x1, y1 = wcs.world_to_pixel(c1)
        x2, y2 = wcs.world_to_pixel(c2)
        sign = (x2 - x1) * (yy - y1) - (y2 - y1) * (xx - x1)

        # Combine them and create an overlay
        circle_mask = circle_mask & (sign > 0)
        # overlay = np.zeros((*side_mask.shape, 4))
        # overlay[..., 3] = 0.0  # fully transparent everywhere
        # overlay[side_mask] = [0, 0, 0, 1]  # black & opaque where mask is True

    return circle_mask


def move_files(file_list, src_root_path, dest_path):

    # Clear the folder if it already exists:
    if dest_path.exists():
        shutil.rmtree(dest_path)

    # Recreate the folder
    dest_path.mkdir(parents=True, exist_ok=True)

    for src in file_list:
        src_path = src_root_path / src
        if src_path.exists():
            print(src_path, '->', dest_path / src_path.name)
            shutil.copy(src_path, dest_path / src_path.name)
        else:
            print(f"-------- File not found: {src_path}")

    return


def add_opacity_profile(spec, opacity_pname, kernel=20, pxs=1.570539778139733, voigtfit_pname=None):

    # Check if true
    if not opacity_pname.is_file():
        print(f'- WARNING: No opacity lines frame at {opacity_pname}')
        spec.plot.show()

    else:
        norm_flux_profile = absorption_spectrum(opacity_pname, spec, kernel=kernel, pxs=pxs)

        # # Container for the opacity
        # tau_arr = np.zeros(spec.wave.data.size)
        #
        # # Loop though the lines and add up the opacities
        # lines_df = lime.load_frame(opacity_pname)
        # for row in lines_df.itertuples():
        #     print(row.Index, f'idx={int(row.idx_comp)}', row.wavelength, row.osc_str, row.trans_prob, row.logN, row.b, row.z_line)
        #     tau_arr += optical_depth_profile(spec.wave.data,
        #                                      row.wavelength,
        #                                      row.osc_str,
        #                                      row.trans_prob,
        #                                      np.power(10, row.logN),
        #                                      row.b,
        #                                      row.z_line)
        #
        # norm_flux_profile = profile_normflux(tau_arr, spec.wave.data, kernel=kernel, pxs=pxs)

        # spec.plot.ax.step(spec.wave, norm_flux_profile, linestyle='--', color=lime.theme.colors['cont'], where='mid',
        #                   label='LiMe')

        if voigtfit_pname.is_file():
            v_wave, v_nflux, v_nerr, v_bestfit, v_mask = np.loadtxt(voigtfit_pname, unpack=True)
            spec.plot.ax.plot(v_wave, v_bestfit, linestyle=':', color='yellow', label='Voigtfit')
            # spec.plot.ax.step(v_wave, v_nflux, linestyle=':', color='yellow', where='mid')

        spec.plot.show()

    return


def on_click(event):

    toolbar = event.canvas.toolbar

    # If a tool (pan/zoom) is active, do nothing
    if toolbar is not None and toolbar.mode != '':
        return

    if event.button == 3 and event.inaxes is not None:
        print(f"x = {event.xdata:.3f}")
        return


class IntervalSelector:
    def __init__(self):
        self.selected_intervals = []

    def on_select(self, xmin, xmax):
        x0, x1 = sorted([xmin, xmax])
        self.selected_intervals.append([x0, x1])
        print(self.selected_intervals)