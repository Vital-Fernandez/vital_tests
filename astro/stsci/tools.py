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



def profile_normflux(tau_arr, wave_arr, kernel, pxs):

    profile_int = np.exp(-tau_arr)

    pad_width = kernel
    padded = np.pad(profile_int, pad_width, mode='edge')  # extend edge values, not zeros

    # Kernel from the LSF
    sigma_instrumental = kernel / 2.35482 / pxs
    LSF = scipy_gaussian(wave_arr.size // 2, sigma_instrumental)
    LSF = LSF / LSF.sum()

    profile_int = fftconvolve(padded, LSF, 'same')
    profile_int = profile_int[pad_width:-pad_width]  #

    # # Perform the convolution
    # fig, ax = plt.subplots()
    # ax.step(wave_arr, profile_int)
    # plt.show()


    return profile_int



def absorption_spectrum(spec, line_list, measurements_dict, kernel=20, pxs=1.570539778139733, min_limit=0.01,
                        plot_fit=False):

    # Container for the opacity
    tau_arr = np.zeros(spec.wave.data.size)

    # Loop though the lines and add up the opacities
    for line in line_list:
        line = Line.from_transition(line)
        for comp, measurements in measurements_dict['components']['HI'].items():
            tau_arr += optical_depth_profile(spec.wave.data,
                                             line.wavelength,
                                             line.atom_data.osc_str,
                                             line.atom_data.trans_prob,
                                             np.power(10, measurements['logN'][0]),
                                             measurements['b_kms'][0],
                                             measurements['redshift'][0])

    norm_spectrum = profile_normflux(tau_arr, spec.wave.data, kernel=kernel, pxs=pxs)
    norm_spectrum[norm_spectrum < min_limit] = 1

    # Plot Normalized spectrum
    if plot_fit:
        spec.plot.spectrum(in_fig=None)
        spec.plot.ax.step(spec.wave, spec.flux/norm_spectrum, linestyle='--', color=lime.theme.colors['cont'], where='mid')
        spec.plot.show()

    return norm_spectrum

def lines_frame_to_lime_lines(df, fit_cfg):

    line_list = []
    print(f'\n- Input lines:')
    for label in df.index:
        line = Line.from_transition(label, data_frame=df, fit_cfg=fit_cfg)
        line_list.append(line)
        msg = f' {line.label} {'' if len(line.list_comps) == 1 else line.list_comps}'

    return line_list


def build_ion_dicts(df, param_hdrs =('b', 'z', 'logN', 'var_z', 'var_b', 'var_N'), default_dict=None):
    """
    Return { ion -> { origin -> { col: value } } } for all ions and
    their observed origins.

    Rules per (ion, origin, column):
      - One unique non-NaN value  → use it.
      - Multiple unique non-NaN values → warn, use the first.
      - All NaN → fall back to DEFAULTS[col].

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns: "ion", "origin", + all keys in DEFAULTS.

    Returns
    -------
    dict  { ion_str: { origin_str: { col: value } } }
    """
    required = {"ion", "origin"} | set()
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"DataFrame is missing columns: {missing}")

    result: dict[str, dict[str, dict]] = {}

    for ion, ion_group in df.groupby("ion", sort=False):
        result[ion] = {}

        for origin, group in ion_group.groupby("origin", sort=False):
            origin_dict = {}

            for col in param_hdrs:
                unique_vals = group[col].dropna().unique()

                if len(unique_vals) == 0:
                    if col.startswith('var'):
                        origin_dict[col] = default_dict[col]
                    else:
                        origin_dict[col] = default_dict[col][ion]

                elif len(unique_vals) == 1:
                    origin_dict[col] = unique_vals[0]

                else:
                    warnings.warn(
                        f"[ion={ion!r}, origin={origin!r}] column {col!r} has "
                        f"{len(unique_vals)} unique values {list(unique_vals)}. "
                        f"Using first: {unique_vals[0]!r}",
                        UserWarning,
                        stacklevel=2,
                    )
                    origin_dict[col] = unique_vals[0]

            result[ion][origin] = origin_dict

    return result


def run_VoigtFit(fpath, spec, lines_df, fit_cfg, conv_dict, resolution=20, obj_redshift=0, output_toml=None,
                 voigt_default_params=None, var_z=True, var_b=True, var_N=True, show_plots=False):

    import VoigtFit as vf

    # Generate spectrum object
    obj_name = fpath.stem
    dataset = vf.DataSet(redshift=obj_redshift, name=obj_name)
    dataset.verbose = True

    dataset.add_data(spec.wave.data, spec.flux.data,
                     res=resolution,
                     err=spec.err_flux.data,
                     mask=~spec.flux.mask,
                     normalized=True,
                     nsub=6)

    # Loop through the lines, add the regions to the data set and check for the fitting transition configuration
    comps_params = {}
    print('\nInput lines:')
    for line in lines_df.index:
        line = lime.Line.from_transition(line, fit_cfg, lines_df)
        msg = f'\n{line}' + f' {line.list_comps if len(line.list_comps) > 1 else ""}'
        print(msg)
        for trans in line.list_comps:
            label_vf = conv_dict[trans.core]

            # Add line with its region
            if label_vf not in dataset.all_lines:
                v_blue = c_KMpS * (line.wavelength - line.mask[2]) / line.wavelength
                v_red = c_KMpS * (line.mask[3] - line.wavelength) / line.wavelength
                print(f'- {trans} = {label_vf} ({-v_blue:0.1f} km/s, {v_red:0.1f} km/s)')
                dataset.add_line(line_tag=label_vf, velspan=(-v_blue, v_red))

            # Kinematic parameters
            kinem = trans.kinem

            # Origin parameters
            orig = 'none' if trans.origin is None else trans.origin
            z_comp = spec.redshift if trans.redshift in [None, 'none'] else trans.redshift

            # Ion based paramters
            ion = trans.classic_notation(just_particle=True)

            # Transition based parameters...
            b = fit_cfg.get(trans, {}).get('b', np.nan)
            logN = fit_cfg.get(trans, {}).get('logN', np.nan)
            var_z_i = fit_cfg.get(trans, {}).get('var_z')
            var_b_i = fit_cfg.get(trans, {}).get('var_b')
            var_N_i = fit_cfg.get(trans, {}).get('var_N')
            comp_label = f'{ion}_{kinem}_{orig}_{z_comp}'
            comps_params[trans] = dict(vf_label=label_vf,
                                       lime_label=comp_label,
                                       line=line,
                                       origin=orig,
                                       kinem=kinem,
                                       ion=ion,
                                       b=b,
                                       z=np.nan if z_comp is None else z_comp,
                                       logN=logN,
                                       var_z=var_z_i,
                                       var_b=var_b_i,
                                       var_N=var_N_i)


    # print(df_lines.to_string())

    # Prepare default values
    voigt_default_params = {**voigt_default_params}
    voigt_default_params['var_z'] = var_z
    voigt_default_params['var_b'] = var_b
    voigt_default_params['var_N'] = var_N

    # Determine components and resolve conflicts
    df_lines = pd.DataFrame(index=comps_params.keys(), data=comps_params.values())
    df_components = build_ion_dicts(df_lines, default_dict=voigt_default_params)

    # Check for conflicts and components
    print(f'- Components:')
    for ion in df_components.keys():
        print(f'-- {ion}: ')
        for origin, orig_params in df_components[ion].items():
            print(f'--- {origin}: ',", ".join(f"{k}: {v}" for k, v in orig_params.items()))
            dataset.add_component(ion=ion, **orig_params)

    print('\nInput Voigtfit lines')
    print(dataset.lines)

    print('\nInput Voigtfit components')
    print(dataset.components)

    # Normalize and mask if necessary
    dataset.prepare_dataset(norm=False, mask=False)

    # Fit the data
    popt, chi2 = dataset.fit(verbose=True, factor=10.)
    dataset.save_fit_regions(filename=str(fpath.parent/f'{obj_name}_best_fit'))

    # Save to toml
    param_dict = extract_voigtfit_params(dataset, obj_name)
    lime.save_cfg(output_toml, param_dict, section_name=f'{obj_name}_results', clear_section=True)
    strip_all_quotes(output_toml)

    # Save the results
    if show_plots:
        dataset.plot_fit(individual=True, xunit='wl')
        plt.show(block=True)
    else:
        dataset.plot_fit(filename=str(fpath.parent / f'{obj_name}_voigtfit_plot.pdf'), individual=False, xunit='wl')

    # Generate fit plot
    # plot_profile_components(str(fpath.parent/f'{obj_name}_comps_fitting'), obj_name, objSpec, dataset)

    # -- Print total column densities
    dataset.print_total()
    dataset.save_parameters(str(fpath.parent/f'{obj_name}_voigtfit_parameters'))
    # dataset.save_cont_parameters_to_file(str(fpath.parent/f'{obj_name}_voigtfit_cont'))

    return


def run_VoigtFit_orig(fpath, objSpec, line_list, comps_dict, conv_dict, resolution=20, obj_redshift=0, output_toml=None):

    import VoigtFit as vf

    # Generate spectrum object
    obj_name = fpath.stem
    dataset = vf.DataSet(redshift=obj_redshift, name=obj_name)
    dataset.verbose = True

    dataset.add_data(objSpec.wave.data, objSpec.flux.data,
                     res=resolution,
                     err=objSpec.err_flux.data,
                     mask=~objSpec.flux.mask,
                     normalized=True,
                     nsub=6)

    # Declare the lines
    print('\nInput lines')
    for i, line in enumerate(line_list):
        label_vf = conv_dict[line.core]
        v_blue = c_KMpS * (line.wavelength - line.mask[2]) / line.wavelength
        v_red = c_KMpS * (line.mask[3] - line.wavelength) / line.wavelength
        msg = (f' {line.label} = {label_vf}  {'' if len(line.list_comps) == 1 else line.list_comps} '
               f'({-v_blue:0.1f} km/s, {v_red:0.1f} km/s)')
        print(msg)
        dataset.add_line(line_tag=label_vf, velspan=(-v_blue, v_red))

    # Define kinematic components from reference line
    if 'kinem_comps' in comps_dict:
        for line_label in comps_dict['kinem_comps'].keys():
            if line_label in line_list:
                line_ref = line_list[line_list.index(line_label)]
                for j, trans in enumerate(line_ref.list_comps):
                    kin_params = {}
                    kin_params['ion'] = trans.particle.classic_notation
                    kin_params['z'] = objSpec.redshift if trans.redshift is None else 0
                    for i, (param, value_list) in enumerate(comps_dict['kinem_comps'][line_ref.label].items()):
                        kin_params[param] = value_list[j]
                    dataset.add_component(**kin_params)
                    print(line_label, kin_params)

    # For all lines
    else:
        for line in line_list:
            for j, trans in enumerate(line.list_comps):
                # ion = trans.particle.classic_notation
                ion = trans.classic_notation(just_particle=True)

                # Check the ion data has not been introduced already
                if len(dataset.components[ion]) == 0:
                    z = objSpec.redshift if trans.redshift is None else 0
                    b = comps_dict.get('b', {}).get(ion, DEFAULT_b[ion])
                    logN = comps_dict.get('logN', {}).get(ion, DEFAULT_logN[ion])
                    var_z = comps_dict.get('var_z', {}).get(ion, True)
                    var_b = comps_dict.get('var_b', {}).get(ion, True)
                    var_N = comps_dict.get('var_N', {}).get(ion, True)
                    kin_params = {'ion':ion, 'z':z, 'b':b, 'logN':logN, 'var_z':var_z, 'var_b':var_b, 'var_N':var_N}
                    dataset.add_component(**kin_params)
                    print(line, kin_params)

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

    print('\nInput Voigtfit lines')
    print(dataset.lines)

    print('\nInput Voigtfit components')
    print(dataset.components)


    # Normalize and mask if necessary
    dataset.prepare_dataset(norm=False, mask=False)

    # Fit the data
    popt, chi2 = dataset.fit(verbose=True, factor=10.)
    dataset.save_fit_regions(filename=str(fpath.parent/f'{obj_name}_best_fit'))

    # Save to toml
    param_dict = extract_voigtfit_params(dataset, obj_name)
    lime.save_cfg(output_toml, param_dict, section_name=f'{obj_name}_results', clear_section=True)
    strip_all_quotes(output_toml)

    # Save the results
    dataset.plot_fit(filename=str(fpath.parent/f'{obj_name}_voigtfit_plot.pdf'), individual=False, xunit='wl')
    # dataset.plot_fit(individual=True, xunit='wl')
    # plt.show(block=True)

    # Generate fit plot
    # plot_profile_components(str(fpath.parent/f'{obj_name}_comps_fitting'), obj_name, objSpec, dataset)

    # -- Print total column densities
    dataset.print_total()
    dataset.save_parameters(str(fpath.parent/f'{obj_name}_voigtfit_parameters'))
    # dataset.save_cont_parameters_to_file(str(fpath.parent/f'{obj_name}_voigtfit_cont'))

    return


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

    norm_profile = profile_normflux(tau_arr, spec.wave.data, kernel=20, pxs=1.570539778139733)
    spec.plot.ax.plot(spec.wave.data, norm_profile, linestyle='--', linewidth=2, label=label)

    spec.plot.ax.legend(loc=9)
    plt.tight_layout()
    plt.savefig(output_fpath)
    np.savetxt(Path(output_fpath).with_suffix('.txt'), np.column_stack([spec.wave.data, norm_profile]))


    return


def extract_voigtfit_params(dataset, prefix, velocity=False):
    """
    Convert best-fit parameters to a single dictionary section
    named after the variable prefix.
    """
    z_sys = dataset.redshift

    # Everything goes into this one dictionary
    combined_data = {
        "z_sys": float(z_sys),
        "chi2": float(dataset.chi2),
        "n_free": int(dataset.minimizer.result.nfree)
    }

    # Add the components to the same dictionary
    for ion in sorted(dataset.components.keys()):
        for i in range(len(dataset.components[ion])):
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

            # Dot notation keys added directly to the main dict
            combined_data[f"components.{ion}.{i}.{param_name}"] = [float(val), float(err)]
            combined_data[f"components.{ion}.{i}.b_kms"] = [float(b.value), float(b.stderr) if b.stderr else 0.0]
            combined_data[f"components.{ion}.{i}.logN"] = [float(logN.value), float(logN.stderr) if logN.stderr else 0.0]

    # Return as a nested dict so TOML creates one header: [prefix]
    return combined_data

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