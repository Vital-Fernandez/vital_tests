import lime
import numpy as np
import pandas as pd
import configparser
import pyneb as pn
import sys

from pathlib import Path
from astropy.io import fits
from mpdaf.obj import Cube
from astropy.table import Table
from lime.tools import label_decomposition
from lmfit.models import LinearModel
from lime.io import progress_bar
from collections import Sequence
from uncertainties import umath, unumpy, ufloat
from lime.plots import STANDARD_PLOT
from matplotlib import pyplot as plt, cm, colors, rc_context
from lime.io import save_cfg
from astropy.wcs import WCS

# State target lines and parameters
target_lines = ['H1_4861A', 'H1_4861A_w1', 'H1_6563A',  'H1_6563A_w1', 'H1_6563A_w2',
                'H1_8750A', 'H1_8863A', 'H1_9015A', 'H1_9229A',
                'O3_4959A', 'O3_5007A', 'O3_5007A_w1',
                'S2_6716A', 'S2_6731A',
                'S3_6312A', 'S3_9069A',
                'He1_5876A', 'He1_6678A',  'He1_4922A', 'He1_5016A', 'He1_6678A', 'He1_7065A', 'He1_7281A', 'He1_8446A']

# State target lines and parameters
abs_target_lines = ['H1_4861A', 'He1_4922A', 'He1_5016A',
                    'He1_5876A', 'H1_6563A', 'He1_6678A',
                    'H1_8750A', 'H1_8863A', 'H1_9015A', 'H1_9229A']

# Generate the parameter maps data
param_images = {'intg_flux': target_lines,
                'intg_err': target_lines,
                'gauss_flux': target_lines,
                'gauss_err': target_lines,
                'v_r': target_lines,
                'v_r_err': target_lines,
                'z_line': target_lines,
                'center': target_lines,
                'center_err': target_lines,
                'sigma_vel': target_lines,
                'sigma_vel_err': target_lines,
                'cont': target_lines}

param_images_abs = {'intg_flux': abs_target_lines,
                    'intg_err': abs_target_lines,
                    'gauss_flux': abs_target_lines,
                    'gauss_err': abs_target_lines,
                    'v_r': abs_target_lines,
                    'v_r_err': abs_target_lines,
                    'center': abs_target_lines,
                    'center_err': abs_target_lines,
                    'sigma_vel': abs_target_lines,
                    'sigma_vel_err': abs_target_lines,
                    'cont': abs_target_lines}

label_Conver = {'H1_6563A': 'Halpha',
                'H1_4861A': 'Hbeta',
                'H1_9229A': 'HPas9',
                'H1_9015A': 'HPas10',
                'H1_8863A': 'HPas11',
                'H1_8750A': 'HPas12'}

lineAreas = {'H1_6563A': (6558.0, 6568.0),
             'S3_6312A': (6308.15, 6317.25),
             'O3_5007A': (5002.0, 5013.0),
             'S2_6717A': (6717.0, 6734.0)}

dinamicLines = {'H1_6563A': r'$H\alpha_{Narrow}$',
              # 'H1_6563A_w1': r'$H\alpha_{Broad\,1}$',
              # 'H1_6563A_w2': r'$H\alpha_{Broad\,2}$',
              'H1_4861A': r'$H\beta_{Narrow}$',
              # 'H1_4861A_w1': r'$H\beta_{Broad\,1}$',
              'O3_4959A': r'$[OIII]5007\AA_{Narrow}$',
              'O3_5007A': r'$[OIII]5007\AA_{Narrow}$',
              # 'O3_5007A_w1': r'$[OIII]5007\AA_{Broad\,1}$',
              'S2_6716A': r'$[SII]6716\AA_{Narrow}$',
              'S2_6731A': r'$[SII]6731\AA_{Narrow}$',
              'S3_9069A': r'$[SIII]9069\AA_{Narrow}$',
              'He1_5876A': r'$HeI\,5876\AA_{Narrow}$'}

latex_Conver = {'H1_6563A': r'H\alpha',
                'H1_4861A': r'H\beta',
                'S2_6731A': r'[SII]6731\AA',
                'H1_9229A': r'H_{Pas,\,9}',
                'H1_9015A': r'H_{Pas,\,10}',
                'H1_8863A': r'H_{Pas,\,11}',
                'H1_8750A': r'H_{Pas,\,12}',
                'v_r': r'$v_{r}\,(km/s)$',
                'sigma_vel': r'$\sigma_{int}\,(km/s)$'}

HII_chemistry_label_conversion = {'logOH': '12+log(O/H)',
                                'logNO': 'log(N/O)',
                                'logU': 'log(U)',
                                'O2_3726A_b': 'OII_3727',
                                'Ne3_3869A': 'NeIII_3868',
                                'O3_4363A': 'OIII_4363',
                                'O3_4959A': 'OIII_4959',
                                'O3_5007A': 'OIII_5007',
                                'N2_6584A': 'NII_6584',
                                'S2_6716A_b': 'SII_6725',
                                'S2_6716A': 'SII_6716',
                                'S2_6731A': 'SII_6731'}


sampling_grid_lines = np.array(['Ar4_4740A',  'O3_4959A',  'O3_5007A', 'N1_5200A',  'Cl3_5518A',  'Cl3_5538A',  'N2_5755A',  'He1_5876A', 'O1_6300A',  'S3_6312A',
                                'N2_6548A',  'H1_6563A',  'N2_6584A',  'He1_6678A',  'S2_6716A',  'S2_6731A',  'He1_7065A',
                                'Ar3_7136A',  'O2_7319A',  'O2_7330A',  'Ar3_7751A',  'S3_9069A',  'S3_9531A'])

HII_chemistry_input_file_columns = ['OIII_4959',  'eOIII_4959',  'OIII_5007',  'eOIII_5007',  'NII_6584',
                                    'eNII_6584',  'SII_6716',  'eSII_6716',  'SII_6731',  'eSII_6731']

latex_labels = {'y_plus': r'$y^{+}$',
             'He1_abund': r'$y^{+}$',
             'He2_abund': r'$y^{++}$',
             'Te': r'$T_{e}$',
             'T_low': r'$T_{low}(K)$',
             'T_LOW': r'$T_{low}(K)$',
             'T_high': r'$T_{high}(K)$',
             'T_HIGH': r'$T_{high}(K)$',
             'T_He': r'$T_{He}$',
             'n_e': r'$n_{e}(cm^{-3})$',
             'cHbeta': r'$c(H\beta)$',
             'tau': r'$\tau$',
             'xi': r'$\xi$',
             'ChiSq': r'$\chi^{2}$',
             'ChiSq_Recomb': r'$\chi^{2}_{Recomb}$',
             'ChiSq_Metals': r'$\chi^{2}_{Metals}$',
             'ChiSq_O': r'$\chi^{2}_{O}$',
             'ChiSq_S': r'$\chi^{2}_{S}$',
             'S2_abund': r'$S^{+}$',
             'He1r': r'$y^{+}$',
             'He2r': r'$y^{2+}$',
             'He1': r'$y^{+}$',
             'He2': r'$y^{2+}$',
             'log(He1r)': r'$log(y^{+})$',
             'log(He2r)': r'$log(y^{2+})$',
             'OH': r'$\frac{O}{H}$',
             'OH_err': r'$O/H\,err$',
             'S3_abund': r'$S^{2+}$',
             'O2_abund': r'$O^{+}$',
             'O3_abund': r'$O^{2+}$',
             'S3_abund': r'$S^{2+}$',
             'O2_abund': r'$O^{+}$',
             'O3_abund': r'$O^{2+}$',
             'N2_abund': r'$N^{+}$',
             'Ar3_abund': r'$Ar^{2+}$',
             'Ar4_abund': r'$Ar^{3+}$',
             'S2': r'$\frac{S^{+}}{H^{+}}$',
             'S3': r'$\frac{S^{2+}}{H^{+}}$',
             'S4': r'$\frac{S^{3+}}{H^{+}}$',
             'O2': r'$\frac{O^{+}}{H^{+}}$',
             'O3': r'$\frac{O^{2+}}{H^{+}}$',
             'Ni3': r'$\frac{Ni^{2+}}{H^{+}}$',
             'NI3': r'$\frac{Ni^{2+}}{H^{+}}$',
             'Cl3': r'$\frac{Cl^{2+}}{H^{+}}$',
             'CL3': r'$\frac{Cl^{2+}}{H^{+}}$',
             'Ne3': r'$\frac{Ne^{2+}}{H^{+}}$',
             'NE3': r'$\frac{Ne^{2+}}{H^{+}}$',
             'Fe3': r'$\frac{Fe^{2+}}{H^{+}}$',
             'FE3': r'$\frac{Fe^{2+}}{H^{+}}$',
             'N2': r'$\frac{N^{+}}{H^{+}}$',
             'Ar3': r'$\frac{Ar^{2+}}{H^{+}}$',
             'AR3': r'$\frac{Ar^{2+}}{H^{+}}$',
             'Ar4': r'$\frac{Ar^{3+}}{H^{+}}$',
             'AR4': r'$\frac{Ar^{3+}}{H^{+}}$',
             'Cl4': r'$\frac{Cl^{3+}}{H^{+}}$',
             'CL4': r'$\frac{Cl^{3+}}{H^{+}}$',
             'Ar_abund': r'$\frac{ArI}{HI}$',
             'He_abund': r'$\frac{HeI}{HI}$',
             'O_abund': r'$\frac{OI}{HI}$',
             'N_abund': r'$\frac{NI}{HI}$',
             'S_abund': r'$\frac{SI}{HI}$',
             'Ymass_O': r'$Y_{O}$',
             'Ymass_S': r'$Y_{S}$',
             'Ar': r'$\frac{Ar}{H}$',
             'He': r'$\frac{He}{H}$',
             'O': r'$\frac{O}{H}$',
             'N': r'$\frac{N}{H}$',
             'S': r'$\frac{S}{H}$',
             'Ymass_O': r'$Y_{O}$',
             'Ymass_S': r'$Y_{S}$',
             'NO': r'$\frac{N}{O}$',
             'calcFluxes_Op': 'Line fluxes',
             'z_star': r'$z_{\star}$',
             'sigma_star': r'$\sigma_{\star}$',
             'Av_star': r'$Av_{\star}$',
             'chiSq_ssp': r'$\chi^{2}_{SSP}$',
             'x': r'x interpolator$',
             'ICF_SIV': r'$ICF\left(S^{3+}\right)$',
             'logU': r'$log(U)$',
             'logOH': r'$\frac{O}{H}$',
             'logNO': r'$log(\frac{N}{O})$',
             'Teff': r'$T_{eff}$',
             'TEFF': r'$T_{eff}$',
             'X_i+': r'$X^{i+}$',
             'NH':   r'$\frac{N}{H}$',
             'ArH': r'$\frac{Ar}{H}$',
             'SH': r'$\frac{S}{H}$',
             'ICF_S4': r'$ICF(S^{3+})$',
             'SO': r'$\frac{S}{O}$',

            'O2_O3': r'$\frac{O^{+}}{O^{2+}}$',
            'S2_S3': r'$\frac{S^{+}}{S^{2+}}$',
             'log(X_i+)': r'$12+log\left(X^{i+}\right)$',
             'redNoise': r'$\Delta(cH\beta)$',
             'nSII_cHbeta': r'$\frac{n_{e,\,[SII]}}{c(H\beta)}$',
             'Y_O': r'$Y_{\frac{O}{H}}$',
             'Y_S': r'$Y_{\frac{S}{H}}$',
             'eta':r'$\eta$'}

signif_figures = {'n_e': 0,
                  'T_low': 0,
                  'T_high': 0,
                  'cHbeta': 2,
                  'Ar4': 2,
                  'Ar3': 2,
                  'O2': 2,
                  'O3': 2,
                  'N2': 2,
                  'He1': 3,
                  'S2': 2,
                  'S3': 2,
                  'OH': 2,
                  'NO': 2,
                  'nSII_cHbeta': 1,
                  'logU': 2,
                  'logOH': 2,
                  'logNO': 2,
                  'O2_O3': 2,
                  'S2_S3': 2,
                  'eta': 2,
                  'OH': 2,
                  'NO': 2,
                  'NH': 2,
                  'ArH': 2,
                  'S4': 2,
                  'SH': 2,
                  'ICF_S4': 2,
                  'SO': 2,
                  'Y_O': 2,
                  'Y_S': 2}

param_units = {'n_e': '$cm^{-3}$',
               'T_low': '$K$',
               'cHbeta': '',
               'Ar4':    '',
               'Ar3':    '',
               'O2':     '',
               'O3':     '',
               'N2':     '',
               'He1':    '',
               'S2':     '',
               'S3':     '',
               'OH':      '',
               'NO':   '',
               'logU': '',
               'logOH': '',
               'logNO': '',
               'nSII_cHbeta': '',
               'O2_O3': '',
               'S2_S3': ''}

convert_dict = {'logOH': 'OH',
                'logNO': 'NO',
                'logU': 'logU'}

technique_convert_label = {'GridSampling': 'Neural model fitting', 'HII-CHI-mistry': r'HII-CHI-mistry', 'direct_method': 'Direct method',
                           'neural_fitting': 'Neural networks'}



def to_natural_abund(abund_array):
    # return np.power(10, abund_array - 12)
    return unumpy.pow(10, abund_array - 12)


def to_log_abund(abund_array):
    # return 12 + np.log10(abund_array)
    return 12 + unumpy.log10(abund_array)


def distribution_parametrisation(param_name, param_array):

    median = np.round(np.nanmedian(param_array), signif_figures[param_name])
    up_lim = np.round(np.nanpercentile(param_array, 84) - np.nanmedian(param_array), signif_figures[param_name])
    low_lim = np.round(np.nanmedian(param_array) - np.nanpercentile(param_array, 16), signif_figures[param_name])
    n_voxels = np.sum(~np.isnan(param_array))

    plot_label = r'{} = ${}^{{{}}}_{{{}}}$ ({} voxels)'.format(latex_labels[param_name], median, up_lim, low_lim, n_voxels)
    log_array = np.array([median, up_lim, low_lim, n_voxels])

    return plot_label, log_array


def plot_parameter_image(plot_db_fits, parameter_list, output_folder, conf_label, tech_label):

    # Image background
    flux6563_image = fits.getdata(plot_db_fits, f'H1_6563A_flux', ver=1)
    halpha_min_level = fits.getval(plot_db_fits, keyword=f'P9050', extname=f'H1_6563A_flux')
    halpha_thresd_level = fits.getval(plot_db_fits, keyword=f'P9250', extname=f'H1_6563A_flux')
    bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)

    # Plot configuration
    defaultConf = STANDARD_PLOT.copy()

    for parameter in parameter_list:

        with fits.open(f'{output_folder}/{conf_label}_{parameter}.fits') as hdul:
            param_image, param_hdr = hdul[parameter].data, hdul[parameter].header

            with rc_context(defaultConf):

                halpha_cmap = cm.gray.copy()
                halpha_cmap.set_under('black')

                fig = plt.figure(figsize=(10, 10))
                ax = fig.add_subplot(projection=WCS(param_hdr), slices=('x', 'y'))

                im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
                im2 = ax.imshow(param_image)
                cbar = fig.colorbar(im2, ax=ax)

                ax.set_title(f'CGCG007−025, {latex_labels[parameter]} \n {technique_convert_label[tech_label]} - {conf_label}')
                ax.set_xlabel(r'RA')
                ax.set_ylabel(r'DEC')
                ax.set_xlim(120, 210)
                ax.set_ylim(110, 220)
                plt.savefig(output_folder/f'{tech_label}_{conf_label}_{parameter}_map')
                # plt.show()
                plt.close(fig)

    return


def compute_parameter_distributions(parameter_list, output_folder, conf_label, mask_file, mask_list, tech_label,
                                    param_conv={}):

    # Distributions containers
    store_dict, err_dict = {}, {}

    # Plot configuration
    defaultConf = STANDARD_PLOT.copy()
    defaultConf['legend.fontsize'] = 16
    defaultConf['figure.figsize'] = (10, 10)

    # Get regions mask:
    spatial_mask_dict = {}
    with fits.open(mask_file) as hdu_masks:
        for mask_name in mask_list:
            mask_data = hdu_masks[mask_name].data.astype(bool)
            spatial_mask_dict[mask_name] = mask_data

    # region_labels = mask_list[:]
    region_idcs_list = list(spatial_mask_dict.values())

    # Loop throught the parameter file images
    for parameter in parameter_list:

        with fits.open(f'{output_folder}/{conf_label}_{parameter}.fits') as hdu_list:

            # Parameter value and error distribution loop
            for map_type in ['', '_err']:

                # Conversion for those times saved in a different format
                fits_hdr_param = param_conv.get(parameter, parameter)

                param_image = hdu_list[f'{fits_hdr_param}{map_type}'].data
                dist_container = []
                label_container = []

                # Plot regions distributions
                with rc_context(defaultConf):

                    for i_mask, mask in enumerate(mask_list):

                        data_dist = param_image[region_idcs_list[i_mask]]
                        label_dist, array_dist = distribution_parametrisation(parameter, data_dist)

                        label_container.append(label_dist)
                        dist_container.append(data_dist)

                        # Store the distribution
                        ref_dist = f'{convert_dict.get(parameter, parameter)}{map_type}_{mask}'
                        store_dict[ref_dist] = array_dist

                    # Make the plot
                    fig, ax = plt.subplots()
                    ax.hist(dist_container, bins=15, label=label_container, stacked=True)
                    ax.set_title(f'CGCG007−025, {latex_labels[parameter]} \n {technique_convert_label[tech_label]} {conf_label} \n regions histogram')
                    ax.set_xlabel(latex_labels[parameter])
                    ax.legend()
                    plt.savefig(output_folder/f'{tech_label}_{conf_label}_{parameter}_regions_histogram{map_type}')
                    plt.close(fig)
                    # plt.show()

                # Plot global distribution
                total_idcs_mask = np.array(region_idcs_list).sum(axis=0).astype(bool)
                ref_global = 'global'

                # Plot total distribution
                with rc_context(defaultConf):

                    data_dist = param_image[total_idcs_mask]
                    label_dist, array_dist = distribution_parametrisation(parameter, param_image)

                    # Store the distribution
                    ref_dist = f'{convert_dict.get(parameter, parameter)}{map_type}_{ref_global}'
                    store_dict[ref_dist] = array_dist

                    # Make the plot
                    fig, ax = plt.subplots()
                    ax.hist(data_dist, bins=15, label=label_dist)
                    ax.set_title(f'CGCG007−025, {latex_labels[parameter]} \n {technique_convert_label[tech_label]} {conf_label} \n all voxels histogram')
                    ax.set_xlabel(latex_labels[parameter])
                    ax.legend()
                    plt.savefig(output_folder / f'{tech_label}_{conf_label}_{parameter}_global_histogram{map_type}')
                    plt.close(fig)

    # Save to a file
    save_cfg('../muse_CGCG007.ini', store_dict, section_name=f'{tech_label}_{conf_label}', clear_section=True)

    return


def deredd_fluxes(obs_flux, obs_err, cHbeta_nom, cHbeta_err, lines_flambda):

    # Generate uncertainty variables to propagate the error
    cHbeta = ufloat(cHbeta_nom, cHbeta_err),
    obsFlux_uarray = unumpy.uarray(obs_flux, obs_err)

    # Compute line intensity
    obsInt_uarray = obsFlux_uarray * unumpy.pow(10, cHbeta * lines_flambda)
    obsInt, obsIntErr = unumpy.nominal_values(obsInt_uarray), unumpy.std_devs(obsInt_uarray)

    return obsInt, obsIntErr


def err_prop_sum(a_err, b_err):
    return np.sqrt(a_err**2 + b_err**2)


def chemical_lines_indexing(input_lines, emis_log, abs_log, chem_cfg, recomb_all=True):

    # Generate new dataframe from
    log = emis_log.copy()
    log.insert(loc=1, column='obsFlux', value=np.nan)
    log.insert(loc=2, column='obsFluxErr', value=np.nan)

    # Use integrated fluxes unless blended lines
    idcs_gaussian = (log.profile_label != 'no')
    log.loc[~idcs_gaussian, 'obsFlux'] = log.loc[~idcs_gaussian, 'intg_flux']
    log.loc[~idcs_gaussian, 'obsFluxErr'] = log.loc[~idcs_gaussian, 'intg_err']
    log.loc[idcs_gaussian, 'obsFlux'] = log.loc[idcs_gaussian, 'gauss_flux']
    log.loc[idcs_gaussian, 'obsFluxErr'] = log.loc[idcs_gaussian, 'gauss_err']

    # Remove the absorption
    if recomb_all == True:
        line_abs = np.array(list(chem_cfg['CGCG007_absorptions'].keys()))
        for i, line in enumerate(line_abs):
            line_abs[i] = line[:line.find('_abs')]
    else:
        line_abs = np.array(['H1_4861A', 'H1_6563A', 'H1_8750A', 'H1_8863A', 'H1_9015A', 'H1_9229A'])

    for line in line_abs:

        if line in log.index:

            emis_flux, emis_err = log.loc[line, 'obsFlux'], log.loc[line, 'obsFluxErr']

            if line in ['H1_6563A', 'H1_4861A']:

                if line in abs_log.index:
                    emis_cont, abs_cont = log.loc[line, 'cont'], abs_log.loc[line, 'cont']
                    norm = emis_cont / abs_cont
                    abs_flux, abs_err = abs_log.loc[line, 'gauss_flux'] * -norm, abs_log.loc[line, 'gauss_err'] * -norm
                else:
                    abs_flux, abs_err = 0.0, 0.0

            else:
                abs_Hbeta = abs_log.loc['H1_4861A', 'gauss_flux']
                emis_cont_Hbeta, abs_cont_Hbeta = log.loc['H1_4861A', 'cont'], abs_log.loc['H1_4861A', 'cont']
                norm = -1 * abs_Hbeta * emis_cont_Hbeta/abs_cont_Hbeta
                abs_flux, abs_err = chem_cfg['CGCG007_absorptions'][f'{line}_abs_array'] * norm


            corr_flux = emis_flux + abs_flux
            corr_err = err_prop_sum(emis_err, abs_err)

            if corr_flux > 0:
                log.loc[line, 'obsFlux'] = corr_flux
                log.loc[line, 'obsFluxErr'] = corr_err

    # Normalized by Hbeta
    flux_Hbeta = log.loc['H1_4861A', 'obsFlux']
    log['obsFlux'] = log['obsFlux'] / flux_Hbeta
    log['obsFluxErr'] = log['obsFluxErr'] / flux_Hbeta

    # N2_6548A missing error
    if 'N2_6548A' in log.index:
        if log.loc['N2_6548A', 'obsFluxErr'] == 0.0:
            log.loc['N2_6548A', 'obsFluxErr'] = log.loc['N2_6584A', 'obsFluxErr']

    # Add up [OII] lines
    if ('O2_7319A' in log.index) and ('O2_7330A' in log.index):
        flux_comb = log.loc['O2_7319A', 'obsFlux'] + log.loc['O2_7330A', 'obsFlux']
        err_comb = np.sqrt(log.loc['O2_7319A', 'obsFluxErr'] ** 2 + log.loc['O2_7330A', 'obsFluxErr'] ** 2)
        log.loc['O2_7319A_b'] = None
        log.loc['O2_7319A_b', ['wavelength', 'obsFlux', 'obsFluxErr', 'observations']] = 7325, flux_comb, err_comb, 'no'
        log.loc['O2_7319A_b', 'ion'] = 'O2'

    # Get indeces of good lines
    idcs_lines = log.index.isin(input_lines) & (log.observations == 'no') & (log.obsFlux > 0) & (log.obsFluxErr > 0)
    lines_remove = log.loc[~idcs_lines].index.values
    log.drop(index=lines_remove, inplace=True)

    return log


def save_log_maps(log_file_address, param_list, output_folder, mask_file_address=None, ext_mask='all',
                    ext_log='_INPUTS', default_spaxel_value=np.nan, output_files_prefix=None, page_hdr={}):

    #TODO add log warning if extension not found
    assert Path(log_file_address).is_file(), f'- ERROR: lines log at {log_file_address} not found'
    assert Path(output_folder).is_dir(), f'- ERROR: Output parameter maps folder {output_folder} not found'
    # Compile the list of voxels to recover the provided masks
    if mask_file_address is not None:

        assert Path(mask_file_address).is_file(), f'- ERROR: mask file at {mask_file_address} not found'

        with fits.open(mask_file_address) as maskHDUs:

            # Get the list of mask extensions
            if ext_mask == 'all':
                if ('PRIMARY' in maskHDUs) and (len(maskHDUs) > 1):
                    mask_list = []
                    for i, HDU in enumerate(maskHDUs):
                        mask_name = HDU.name
                        if mask_name != 'PRIMARY':
                            mask_list.append(mask_name)
                    mask_list = np.array(mask_list)
                else:
                    mask_list = np.array(['PRIMARY'])
            else:
                mask_list = np.array(ext_mask, ndmin=1)

            # Combine all the mask voxels into one
            for i, mask_name in enumerate(mask_list):
                if i == 0:
                    mask_array = maskHDUs[mask_name].data
                    image_shape = mask_array.shape
                else:
                    assert image_shape == maskHDUs[
                        mask_name].data.shape, '- ERROR: Input masks do not have the same dimensions'
                    mask_array += maskHDUs[mask_name].data

            # Convert to boolean
            mask_array = mask_array.astype(bool)

            # List of spaxels in list [(idx_j, idx_i), ...] format
            spaxel_list = np.argwhere(mask_array)

    # No mask file is provided and the user just defines an image size tupple (nY, nX)
    else:
        exit()

    # Generate containers for the data:
    images_dict = {}
    for param in param_list:
        images_dict[f'{param}'] = np.full(image_shape, default_spaxel_value)
        images_dict[f'{param}_err'] = np.full(image_shape, default_spaxel_value)

    # Loop through the spaxels and fill the parameter images
    n_spaxels = spaxel_list.shape[0]
    spaxel_range = np.arange(n_spaxels)

    with fits.open(log_file_address) as logHDUs:

        for i_spaxel in spaxel_range:
            idx_j, idx_i = spaxel_list[i_spaxel]
            spaxel_ref = f'{idx_j}-{idx_i}{ext_log}'

            progress_bar(i_spaxel, n_spaxels, post_text=f'spaxels treated ({n_spaxels})')

            # Confirm log extension exists
            if spaxel_ref in logHDUs:

                # Recover extension data
                log_data, log_header = logHDUs[spaxel_ref].data, logHDUs[spaxel_ref].header

                # Loop through the parameters and the lines:
                for param in param_list:
                    if param in log_header:
                        images_dict[f'{param}'][idx_j, idx_i] = log_header[param]
                        images_dict[f'{param}_err'][idx_j, idx_i] = log_header[f'{param}_err']

    # New line after the rustic progress bar
    print()

    # Save the parameter maps as individual fits files with one line per page
    output_files_prefix = '' if output_files_prefix is None else output_files_prefix
    for param in param_list:

        # Primary header
        paramHDUs = fits.HDUList()
        paramHDUs.append(fits.PrimaryHDU())

        # ImageHDU for the parameter maps
        hdr = fits.Header({'PARAM': param})
        hdr.update(page_hdr)
        data = images_dict[f'{param}']
        paramHDUs.append(fits.ImageHDU(name=param, data=data, header=hdr, ver=1))

        # ImageHDU for the parameter error maps
        hdr = fits.Header({'PARAMERR': param})
        hdr.update(page_hdr)
        data_err = images_dict[f'{param}_err']
        paramHDUs.append(fits.ImageHDU(name=f'{param}_err', data=data_err, header=hdr, ver=1))

        # Write to new file
        output_file = Path(output_folder) / f'{output_files_prefix}{param}.fits'
        paramHDUs.writeto(output_file, overwrite=True, output_verify='fix')

    return


def total_abundances_calculation(param_list, output_folder, mask_file_address, regions_list='all', ref_fits='', image_size=None,
                                 header={}):

    # Containers for the parameters
    store_dict, voxels_dict = {}, {}

    # Get regions mask:
    spatial_mask_dict = {}
    with fits.open(mask_file_address) as hdu_masks:
        hdulist = regions_list if regions_list != all else list(hdu_masks.keys())
        for mask_name in hdulist:
            mask_data = hdu_masks[mask_name].data.astype(bool)
            spatial_mask_dict[mask_name] = mask_data
        if image_size is None:
            image_size = hdu_masks[regions_list[0]].data.shape

    region_idcs_list = list(spatial_mask_dict.values())
    total_idcs_mask = np.array(region_idcs_list).sum(axis=0).astype(bool)

    for param in param_list:

        with fits.open(f'{output_folder}/{ref_fits}{param}.fits') as hdu_list:

            image_data, err_data = hdu_list[param].data, hdu_list[f'{param}_err'].data

            array_data = image_data[total_idcs_mask]
            array_err = err_data[total_idcs_mask]

            voxels_dict[param] = unumpy.uarray(array_data, array_err)

    T_low = voxels_dict['T_low']
    T_high = 0.8403 * T_low + 2689

    # Oxygen abundance
    O2, O3 = to_natural_abund(voxels_dict['O2']), to_natural_abund(voxels_dict['O3'])
    OH = O2 + O3

    # Nirogen abundance
    N2 = to_natural_abund(voxels_dict['N2'])
    NO = N2/O2
    NH = NO * OH

    # Argon abundance
    Ar3, Ar4 = to_natural_abund(voxels_dict['Ar3']), to_natural_abund(voxels_dict['Ar4'])
    ArH = Ar3 + Ar4

    # Sulfur abundance
    S2, S3 = to_natural_abund(voxels_dict['S2']), to_natural_abund(voxels_dict['S3'])
    # m_conf, n_conf = np.random.normal(1.162, 0.006, Ar3.size), np.random.normal(0.05, 0.01, Ar3.size)
    m_coef, n_coef = ufloat(1.162, 0.006), ufloat(0.05, 0.01)
    exp_value = (unumpy.log10(Ar3/Ar4) - n_coef) / m_coef
    S3S4 = unumpy.pow(10, exp_value)
    S4 = S3/S3S4
    SH = S2 + S3 + S4
    ICF_S4 = SH/(S2 + S3)
    SO = SH/OH

    #Extra params
    O2_O3 = O2/O3
    S2_S3 = S2/S3
    eta = O2_O3/S2_S3

    # Save mean values to log
    He1 = voxels_dict['He1']
    SO_dors = unumpy.pow(10, -1.78), unumpy.pow(10, -0.02)
    # SO_dist = np.random.normal(SO_dors[0], SO_dors[1], He1.size)
    SO_coeff = ufloat(SO_dors[0], SO_dors[1])
    Y_O = (4 * He1 * (1 - 20 * OH)) / (1 + 4*He1)
    Y_S = (4 * He1 * (1 - 20 * SO_coeff * SH)) / (1 + 4*He1)

    # Store to a dictiory
    param_chain_dict = dict(T_high=T_high,
                            OH=to_log_abund(OH),
                            NO=unumpy.log10(NO),
                            NH=to_log_abund(NH),
                            ArH=to_log_abund(ArH),
                            S4=to_log_abund(S4),
                            SH=to_log_abund(SH),
                            ICF_S4=ICF_S4,
                            SO=unumpy.log10(SO),
                            Y_O=Y_O,
                            Y_S=Y_S,
                            S2_S3=S2_S3,
                            O2_O3=O2_O3,
                            eta=eta)

    # Saving the total abundances as dictionaries
    for param, voxels_array in param_chain_dict.items():

        value_array, error_array = unumpy.nominal_values(voxels_array), unumpy.std_devs(voxels_array)

        param_image = np.full(image_size, np.nan)
        err_image = np.full(image_size, np.nan)
        param_image[total_idcs_mask] = value_array
        err_image[total_idcs_mask] = error_array

        # Primary header
        paramHDUs = fits.HDUList()
        paramHDUs.append(fits.PrimaryHDU())

        # ImageHDU for the parameter maps
        hdr = fits.Header({'PARAM': param})
        hdr.update(header)
        paramHDUs.append(fits.ImageHDU(name=param, data=param_image, header=hdr, ver=1))
        paramHDUs.append(fits.ImageHDU(name=f'{param}_err', data=err_image, header=hdr, ver=1))

        # Write to new file
        output_file = Path(output_folder)/f'{ref_fits}{param}.fits'
        paramHDUs.writeto(output_file, overwrite=True, output_verify='fix')

    return


def muse_grid_sampling_fluxes(lines_array, log, ext_coef, R_V, red_law):

    log.insert(loc=1, column='obsFlux', value=np.nan)
    log.insert(loc=2, column='obsFluxErr', value=np.nan)
    log.insert(loc=3, column='obsInt', value=np.nan)
    log.insert(loc=4, column='obsIntErr', value=np.nan)

    idcs_gaussian = (log.blended_label != 'None') & (~log.index.str.contains('_m'))
    log.loc[~idcs_gaussian, 'obsFlux'] = log.loc[~idcs_gaussian, 'intg_flux']
    log.loc[~idcs_gaussian, 'obsFluxErr'] = log.loc[~idcs_gaussian, 'intg_err']
    log.loc[idcs_gaussian, 'obsFlux'] = log.loc[idcs_gaussian, 'gauss_flux']
    log.loc[idcs_gaussian, 'obsFluxErr'] = log.loc[idcs_gaussian, 'gauss_err']

    # Normalize by Hbeta the lines log for the fitting
    flux_Hbeta = log.loc['H1_4861A', 'obsFlux']
    log['obsFlux'] = log['obsFlux'] / flux_Hbeta
    log['obsFluxErr'] = log['obsFluxErr'] / flux_Hbeta
    ion_array, wave_array, latex_array = lime.label_decomposition(log.index.values)

    # Correct fluxes from extinction
    redCorr = pn.RedCorr(R_V=R_V, law=red_law, cHbeta=ext_coef)
    corr = redCorr.getCorrHb(wave_array)
    log['obsInt'] = log['obsFlux'] * corr
    log['obsIntErr'] = log['obsFluxErr'] * corr

    # Slice to the target lines
    idcs_lines = log.index.isin(lines_array)
    output_log = log.loc[idcs_lines]

    return output_log


def voxel_security_check(linesDF):

    check = False

    if 'S3_6312A' in linesDF.index:
        check = True

    return check


def store_frame_to_fits(fits_address, fits_hdu, ext_name):

    if fits_address.is_file():
        try:
            fits.update(fits_address, data=fits_hdu.data, header=fits_hdu.header, extname=ext_name, verify=True)
        except KeyError:
            fits.append(fits_address, data=fits_hdu.data, header=fits_hdu.header, extname=ext_name)
    else:
        fits_hdu.writeto(fits_address, overwrite=True, output_verify='fix')

    return


def import_muse_fits(file_address):
    cube = Cube(filename=str(file_address))
    header = cube.data_header

    cube.wave.info()
    dw = header['CD3_3']
    w_min = header['CRVAL3']
    nPixels = header['NAXIS3']
    w_max = w_min + dw * nPixels
    wave = np.linspace(w_min, w_max, nPixels, endpoint=False)

    return wave, cube, header


def import_fado_cube(file_address, ext=0):

    with fits.open(file_address) as hdu_list:

        data = hdu_list[ext].data
        hdr = hdu_list[ext].header

        dw = hdr['CDELT3']
        w_min = hdr['CRVAL3']
        nPixels = hdr['NAXIS3']
        w_max = w_min + dw * nPixels
        wave = np.linspace(w_min, w_max, nPixels, endpoint=False)

    return wave, data, hdr


def read_lines_fits(files_dict, night_list):

    lineLabel_conver = {'H19_3686': 'H1_3686A',
                        'H18_3691': 'H1_3691A',
                        'H17_3697': 'H1_3697A',
                        'H16_3703': 'H1_3703A',
                        'H15_3711': 'H1_3711A',
                        'H14_3721': 'H1_3721A',
                        '[OII]_3726_narrow': 'O2_3726A_N1',
                        '[OII]_3726_medium': 'O2_3726A_M1',
                        '[OII]_3728_narrow': 'O2_3729A_N1',
                        '[OII]_3728_medium': 'O2_3729A_M1',
                        'H13_3734': 'H1_3734A',
                        'H12_3750': 'H1_3750A',
                        'H11_3770': 'H1_3770A',
                        'H10_3797': 'H1_3797A',
                        'H9_3835': 'H1_3835A',
                        '[NeIII]_3868_narrow': 'Ne3_3868A_N1',
                        '[NeIII]_3868_medium': 'Ne3_3868A_M1',
                        'HeI_3888': 'He1_3888A',
                        'H8_3889_narrow': 'H1_3889A_N1',
                        'H8_3889_medium': 'H1_3889A_M1',
                        '[NeIII]_3967_narrow': 'Ne3_3967A_N1',
                        '[NeIII]_3967_medium': 'Ne3_3967A_M1',
                        'He_3970_narrow': 'He1_3970A_N1',
                        'He_3970_medium': 'He1_3970A_M1',
                        'HeI_4026': 'He1_4026A',
                        '[SII]_4068': 'S2_4069A',
                        '[SII]_4076': 'S2_4076A',
                        'Hd_4101_narrow': 'H1_4101A_N1',
                        'Hd_4101_medium': 'H1_4101A_M1',
                        'HeI_4143': 'He1_4143A',
                        'Hc_4340_narrow': 'H1_4340A_N1',
                        'Hc_4340_medium': 'H1_4340A_M1',
                        '[OIII]_4363_narrow': 'O3_4363A_N1',
                        '[OIII]_4363_medium': 'O3_4363A_M1',
                        'HeI_4471_narrow': 'He1_4471A_N1',
                        'HeI_4471_medium': 'He1_4471A_M1',
                        'HeII_4685': 'He2_4686A',
                        '[FeIII]_4658': 'Fe3_4658A',
                        '[ArIV]_4711': 'Ar4_4711A',
                        '[ArIV]_4740_narrow': 'Ar4_4740A_N1',
                        '[ArIV]_4740_medium': 'Ar4_4740A_M1',
                        '[OIII]_4958_a_broad': 'O3_4959A_B1',
                        '[OIII]_4958_c_narrow': 'O3_4959A_N3',
                        '[OIII]_5006_c_narrow': 'O3_5007A_N3',
                        'Hb_4861_a_narrow': 'H1_4861A_N1',
                        'Hb_4861_a_medium': 'H1_4861A_M1',
                        'Hb_4861_b_narrow': 'H1_4861A_N2',
                        '[OIII]_4958_a_narrow': 'O3_4959A_N1',
                        '[OIII]_4958_a_medium': 'O3_4959A_M1',
                        '[OIII]_4958_b_narrow': 'O3_4959A_N2',
                        '[FeIII]_4987': 'Fe3_4987A',
                        '[OIII]_5006_a_narrow': 'O3_5007A_N1',
                        '[OIII]_5006_a_medium': 'O3_5007A_M1',
                        '[OIII]_5006_a_broad': 'O3_5007A_B1',
                        '[OIII]_5006_b_narrow': 'O3_5007A_N2',
                        'HeI_5015': 'He1_5016A',
                        'FeII_5197': 'Fe2_5197A',
                        '[NI]_5200': 'N1_5200A',
                        '[FeIII]_5270': 'Fe3_5270A',
                        '[ClIII]_5517': 'Cl3_5517A',
                        '[ClIII]_5537': 'Cl3_5537A',
                        '[NII]_5754': 'N2_5754A',
                        'HeI_5875_narrow': 'He1_5876A_N1',
                        'HeI_5875_medium': 'He1_5876A_M1',
                        '[OI]_6300': 'O1_6300A',
                        '[SIII]_6312': 'S3_6312A',
                        'SiII_6347': 'Si2_6347A',
                        '[OI]_6363': 'O1_6363A',
                        '[NII]_6548_narrow': 'N2_6548A_N1',
                        '[NII]_6548_medium': 'N2_6548A_M1',
                        'Ha_6562_a_narrow': 'H1_6563A_N1',
                        'Ha_6562_a_medium': 'H1_6563A_M1',
                        'Ha_6562_a_broad': 'H1_6563A_B1',
                        'Ha_6562_b_narrow': 'H1_6563A_N2',
                        '[NII]_6583_narrow': 'N2_6584A_N1',
                        '[NII]_6583_medium': 'N2_6584A_M1',
                        'HeI_6678_narrow': 'He1_6678A_N1',
                        'HeI_6678_medium': 'He1_6678A_M1',
                        '[SII]_6716_narrow': 'S2_6716A_N1',
                        '[SII]_6716_medium': 'S2_6716A_M1',
                        '[SII]_6730_narrow': 'S2_6731A_N1',
                        '[SII]_6730_medium': 'S2_6731A_M1',
                        'OI_7002': 'O1_7002A',
                        'HeI_7065_narrow': 'He1_7065A_N1',
                        'HeI_7065_medium': 'He1_7065A_M1',
                        '[ArIII]_7135_narrow': 'Ar3_7135A_N1',
                        '[ArIII]_7135_medium': 'Ar3_7135A_M1',
                        '[ArIV]_7170': 'Ar4_7170A',
                        '[ArIV]_7237': 'Ar4_7237A',
                        'OI_7254': 'O1_7254A',
                        '[ArIV]_7262': 'Ar4_7262A',
                        'HeI_7281': 'He1_7281A',
                        '[OII]_7319_a_narrow': 'O2_7319A_N1',
                        '[OII]_7319_b_narrow': 'O2_7319A_M1',
                        '[OII]_7330_a_narrow': 'O2_7330A_N1',
                        '[OII]_7330_b_narrow': 'O2_7330A_M1',
                        '[NiII]_7377': 'Ni2_7377A',
                        '[NiII]_7411': 'Ni2_7411A',
                        '[FeII]_7452': 'Fe2_7452A',
                        'NI_7468': 'N1_7468A',
                        '[ArIII]_7751': 'Ar3_7751A',
                        'HeI_7816': 'He1_7861A',
                        'Pa20_8392': 'H1_8329A',
                        'Pa19_8413': 'H1_8413A',
                        'Pa18_8437': 'H1_8437A',
                        'OI_8446': 'O1_8446A',
                        'Pa17_8467': 'H1_8467A',
                        'Pa16_8502': 'H1_8502A',
                        'Pa15_8545': 'H1_8545A',
                        'Pa14_8598': 'H1_8598A',
                        'Pa13_8665': 'H1_8665A',
                        'NI_8680': 'N1_8680A',
                        'NI_8703': 'N1_8703A',
                        'NI_8711': 'N1_8711A',
                        'Pa12_8750': 'H1_8750A',
                        'Pa11_8862': 'H1_8862A',
                        'Pa10_9014': 'H1_9014A',
                        '[SIII]_9068_narrow': 'S3_9069A_N1',
                        '[SIII]_9068_medium': 'S3_9069A_M1'}

    columns_conversion = {'Flux': 'gauss_flux',
                          'FluxError': 'gauss_err',
                          'Center': 'mu',
                          'CenterError': 'mu_err',
                          'Sigma': 'sigma',
                          'SigmaError': 'sigma_err'}

    data_dict = {arm: {i_night: None for i_night in night_list} for arm in files_dict.keys()}

    for color, fits_file in files_dict.items():
        with fits.open(fits_file) as hdul:
            for i in night_list:

                # Fits to dataframe
                linesDF = Table(hdul[i].data).to_pandas()
                linesDF.set_index('LineName', inplace=True)
                linesDF.rename(columns_conversion, axis=1, inplace=True)
                linesDF.rename(index=lineLabel_conver, inplace=True)

                # Store to a dictionary
                data_dict[color][i] = pd.DataFrame.copy(linesDF)

    return data_dict


def compute_cHbeta(line_df, reddening_curve, R_v, temp=10000.0, den=100.0, ref_wave='H1_4861A',
                   compMode='auto'):

    assert ref_wave in line_df.index, f'- ERROR: Reference line {ref_wave} is not in input dataframe index'

    # Create hydrogen recombination atom for emissivities calculation
    H1 = pn.RecAtom('H', 1)

    # Use all the lines from the input data frame
    line_labels = line_df.index.values
    ion_ref, waves_ref, latexLabels_ref = label_decomposition(ref_wave, scalar_output=True)
    ion_array, waves_array, latexLabels_array = label_decomposition(line_labels)

    # Mode 1: Distinguish between single (intg_flux) and  blended (gauss_flux) lines
    if compMode == 'auto':
        Href_flux, Href_err = line_df.loc[ref_wave, 'intg_flux'], line_df.loc[ref_wave, 'intg_err']

        obsFlux, obsErr = np.empty(line_labels.size), np.empty(line_labels.size)
        idcs_intg = (line_df.blended == 'None')

        obsFlux[idcs_intg], obsErr[idcs_intg] = line_df.loc[idcs_intg, ['intg_flux', 'intg_err']].values
        obsFlux[~idcs_intg], obsErr[~idcs_intg] = line_df.loc[~idcs_intg, ['gauss_flux', 'gauss_err']].values


    # Mode 2: Use always the gaussian flux
    elif compMode == 'gauss':
        Href_flux, Href_err = line_df.loc[ref_wave, 'gauss_flux'], line_df.loc[ref_wave, 'gauss_err']
        obsFlux, obsErr = line_df['gauss_flux'].values, line_df['gauss_err'].values


    # Ratio propagating the uncertainty between the lines
    obsFlux_norm = obsFlux / Href_flux
    obsErr_norm = obsFlux_norm * np.sqrt(np.square(obsErr/obsFlux) + np.square(Href_err/Href_flux))

    assert not np.any(np.isnan(obsFlux)) in obsFlux, '- ERROR: nan entry in input fluxes for c(Hbeta) calculation'
    assert not np.any(np.isnan(obsErr)) in obsErr, '- ERROR: nan entry in input uncertainties for c(Hbeta) calculation'

    # Theoretical ratios
    refEmis = H1.getEmissivity(tem=temp, den=den, wave=waves_ref)
    emisIterable = (H1.getEmissivity(tem=temp, den=den, wave=wave) for wave in waves_array)
    linesEmis = np.fromiter(emisIterable, float)
    theoRatios = linesEmis / refEmis

    # Reddening law
    rc = pn.RedCorr(R_V=R_v, law=reddening_curve)
    Xx_ref, Xx = rc.X(waves_ref), rc.X(waves_array)
    f_lines = Xx / Xx_ref - 1
    f_ref = Xx_ref / Xx_ref - 1

    # cHbeta slope fit axes
    x_fred = f_lines - f_ref
    y_flux = np.log10(theoRatios) - np.log10(obsFlux_norm)
    y_err = (obsErr_norm/obsFlux_norm) * (1.0 / np.log(10))

    # Perform fit
    lineModel = LinearModel()
    pars = lineModel.make_params(intercept=y_flux.min(), slope=0)
    output = lineModel.fit(y_flux, pars, x=x_fred, weights=1/np.sqrt(y_err))
    cHbeta, cHbeta_err = output.params['slope'].value, output.params['slope'].stderr
    intercept, intercept_err = output.params['intercept'].value, output.params['intercept'].stderr

    # Store the results
    output_dict = dict(cHbeta=cHbeta,
                       cHbeta_err=cHbeta_err,
                       intercept=intercept,
                       intercept_err=intercept_err,
                       obsRecomb=obsFlux_norm,
                       obsRecombErr=obsErr_norm,
                       y=y_flux, y_err=y_err,
                       x=x_fred,
                       line_labels=latexLabels_array,
                       ref_line=latexLabels_ref)

    return output_dict


def formatStringOutput(value, key, section_label=None, float_format=None, nan_format='nan'):

    # TODO this one should be the default option
    # TODO add more cases for dicts
    # Check None entry
    if value is not None:

        # Check string entry
        if isinstance(value, str):
            formatted_value = value

        else:

            # Case of an array
            scalarVariable = True
            if isinstance(value, (Sequence, np.ndarray)):

                # Confirm is not a single value array
                if len(value) == 1:
                    value = value[0]

                # Case of an array
                else:
                    scalarVariable = False
                    formatted_value = ','.join([str(item) for item in value])

            if scalarVariable:

                # Case single float
                if isinstance(value, str):
                    formatted_value = value
                else:
                    if np.isnan(value):
                        formatted_value = nan_format
                    else:
                        formatted_value = str(value)

    else:
        formatted_value = 'None'

    return formatted_value


def safe_cfg(output_file, param_dict, section_name=None, clear_section=False):

    """
    This function safes the input dictionary into a configuration file. If no section is provided the input dictionary
    overwrites the data

    """

    # Creating a new file (overwritting old if existing)
    if section_name is None:

        # Check all entries are dictionaries
        values_list = [*param_dict.values()]
        section_check = all(isinstance(x, {}) for x in values_list)
        assert section_check, f'ERROR: Dictionary for {output_file} cannot be converted to configuration file. ' \
                              f'Confirm all its values are dictionaries'

        output_cfg = configparser.ConfigParser()
        output_cfg.optionxform = str

        # Loop throught he sections and options to create the files
        for section_name, options_dict in param_dict.items():
            output_cfg.add_section(section_name)
            for option_name, option_value in options_dict.items():
                option_formatted = formatStringOutput(option_value, option_name, section_name)
                output_cfg.set(section_name, option_name, option_formatted)

        # Save to a text format
        with open(output_file, 'w') as f:
            output_cfg.write(f)

    # Updating old file
    else:

        # Confirm file exists
        file_check = Path.is_file(output_file)

        # Load original cfg
        if file_check:
            output_cfg = configparser.ConfigParser()
            output_cfg.optionxform = str
            output_cfg.read(output_file)
        # Create empty cfg
        else:
            output_cfg = configparser.ConfigParser()
            output_cfg.optionxform = str

        # Clear section upon request
        if clear_section:
            if output_cfg.has_section(section_name):
                output_cfg.remove_section(section_name)

        # Add new section if it is not there
        if not output_cfg.has_section(section_name):
            output_cfg.add_section(section_name)

        # Map key values to the expected format and store them
        for option_name, option_value in param_dict.items():
            option_formatted = formatStringOutput(option_value, option_name, section_name)
            output_cfg.set(section_name, option_name, option_formatted)

        # Save to a text file
        with open(output_file, 'w') as f:
            output_cfg.write(f)

    return


#Function for interpolation of grids
def interpolate(grid, z, zmin, zmax, n):

   library_file = '/home/vital/Dropbox/Astrophysics/Tools/HCm_v5.22/Libraries_opt/C17_POPSTAR_1myr.dat'

   #Columns of the library
   n_comments = 0
   with open(library_file, 'r') as file1:
      for line in file1:
         if line[0] == '#':
            n_comments += 1

   auxiliar_labels = np.genfromtxt(library_file, dtype=None, names=True, encoding='ascii', skip_header=n_comments).dtype.names
   ncol = len(auxiliar_labels)
   vec = []

   if z == 2:
      label_z = 'logU'
   if z == 1:
      label_z = 'logNO'
   if z == 0:
      label_z = '12logOH'
   type_list_names = []

   for col in auxiliar_labels:
      inter = 0
      no_inter = 0
      type_list_names.append((col, float))
      for row in range(0,len(grid)):
         if grid[label_z][row] < zmin or grid[label_z][row] > zmax: continue
         if z == 2: x = '12logOH'; y = 'logNO'
         if z == 1: x = '12logOH'; y = 'logU'
         if z == 0: x = 'logNO'; y = 'logU'
         if row == (len(grid)-1):
            vec.append(grid[col][row])
            no_inter = no_inter + 1
         elif grid[x][row] < grid[x][row+1] or grid[y][row] < grid[y][row+1] :
            vec.append(grid[col][row])
            no_inter = no_inter + 1
         else:
            inter = inter + 1
            for index in range(0,n):
               i = grid[col][row]+(index)*(grid[col][row+1]-grid[col][row])/n
               vec.append(i)

   out_aux = np.transpose(np.reshape(vec,(-1,n*inter+no_inter)))
   out = np.zeros(out_aux.shape[0], dtype=type_list_names)

   for col_n in range(0, len(auxiliar_labels)):
      out[auxiliar_labels[col_n]] = out_aux[:, col_n]

   return out


def interpolate_orig(grid, z, zmin, zmax, n, ncol=9):


    # ncol = 9 for HII_CHIM_original, ncol=10 for Teff

    vec = []
    for col in range(ncol):
      inter = 0
      no_inter = 0
      for row in range(0,len(grid)):
         if grid[row,z] < zmin or grid[row,z] > zmax: continue
         if z == 2: x = 0; y = 1
         if z == 1: x = 0; y = 2
         if z == 0: x = 1; y = 2
         if row == (len(grid)-1):
            vec.append(grid[row,col])
            no_inter = no_inter + 1
         elif grid[row,x] < grid[row+1, x] or grid[row, y] < grid[row+1, y]:
            vec.append(grid[row,col])
            no_inter = no_inter + 1
         else:
            inter = inter + 1
            for index in range(0,n):
               i = grid[row, col]+(index)*(grid[row+1, col]-grid[row, col])/n
               vec.append(i)
    out = np.transpose(np.reshape(vec, (-1, n*inter+no_inter)))

    return out

def epm_HII_CHI_mistry_orig(input00, output_file, n, sed, inter, const=1, HCm_folder=None):

    if HCm_folder is None:
        HCm_folder = Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v5.22/')

    # Description of the code
    print('-------------------------------------------------')
    print('This is HII-CHI-mistry v. 5.22')
    print(' See Perez-Montero, E. (2014) and Perez-Montero et al. (2019, 2021) for details')
    print(' Insert the name of your input text file with all or some of the following columns:')
    print('')
    print('[OII] 3727')
    print('[NeIII] 3868')
    print('[OIII] 4363')
    print('[OIII] 4959 or [OIII] 5007')
    print('[NII] 6584')
    print('[SII] 6725 or [SII] 6717 and [SII]6731')
    print('with their corresponding labels and errors in adjacent columns')
    print('-------------------------------------------------')
    print(' ')

    input0 = input00
    input1 = input00

    # POPSTAR MODEL
    if sed == 1:
        file_lib = 'C17_POPSTAR_1myr.dat'
        # Counting comments:
        n_comments = 0
        with open(HCm_folder/'Libraries_opt/' + file_lib, 'r') as file6:
            for line in file6:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. No interpolation'
            print('No interpolation for the POPSTAR models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            print('')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. Interpolation'
            print('Interpolation for the POPSTAR models is going to be used.')
            print('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
            print('')
            res_NO = 0.125

    # BPASS MODEL
    elif sed == 2:
        file_lib = 'C17_BPASS_IMF135_mup300_1myr.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file7:
            for line in file7:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'BPASS a_IMF = 1.35, M_up = 300, age = 1Myr. No interpolation'
            print('No interpolation for the BPASS models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            print('')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'BPASS v.2.1, a_IMF = 1.35, M_up = 300, age = 1Myr. Interpolation'
            print('Interpolation for the BPASS  models is going to be used.')
            print('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
            print('')
            res_NO = 0.125

    # Further questions on the AGN models
    if sed == '3':
        # SLOPE ALPHA
        question = True
        while question:
            if int(sys.version[0]) < 3:
                alpha = input('Choose value for alpha(OX) in the AGN models: [1] -0.8 [2] -1.2: ')
            else:
                alpha = input('Choose value for alpha(OX) in the AGN models: [1] -0.8 [2] -1.2:  ')
            if alpha == '1' or alpha == '2': question = False
            alpha = int(alpha)
        print('')
        # Fraction of free electrons (stopping criteria in the models)
        question = True
        while question:
            if int(sys.version[0]) < 3:
                efrac = input(
                    'Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons: ')
            else:
                efrac = input(
                    'Choose stop criterion in the AGN models: [1] 2% free electrons [2] 98% free electrons:  ')
            if efrac == '1' or efrac == '2': question = False
            efrac = int(efrac)
        print('')

    # AGN MODEL FOR alpha_OX = -0.8, efrac = 2% and logU in [-4.0, -0.5]
    elif sed == 3 and alpha == 1 and efrac == 1:
        file_lib = 'C17_AGN_alpha08_efrac02_CNfix.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file8:
            for line in file8:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 2%. No interpolation. No constraint in ionization parameter.'
            print('No interpolation for the AGN a(ox) = -0.8 with 2% free electrons models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O. ')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 2% interpolated. No constraint in ionization parameter.'
            print('Interpolation for the AGN a(ox) = -0.8 with 2% free electrons models is going to be used.')
            print('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
            res_NO = 0.125

    # AGN MODEL FOR alpha_OX = -0.8, efrac = 98% and logU in [-4.0, -0.5]
    elif sed == 3 and alpha == 1 and efrac == 2:
        file_lib = 'C17_AGN_alpha08_efrac98_CNfix.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file9:
            for line in file9:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 98%. No interpolation. No constraint in ionization parameter.'
            print('No interpolation for the AGN a(ox) = -0.8 with 98% free electrons models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 98% interpolated. No constraint in ionization parameter.'
            print('Interpolation for the AGN a(ox) = -0.8 with 98% free electrons models is going to be used.')
            print('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
            res_NO = 0.125

    # AGN MODEL FOR alpha_OX = -1.2, efrac = 2% and logU in [-4.0, -0.5]
    elif sed == 3 and alpha == 2 and efrac == 1:
        file_lib = 'C17_AGN_alpha12_efrac02_CNfix.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file10:
            for line in file10:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 2%. No interpolation. No constraint in ionization parameter.'
            print('No interpolation for the AGN a(ox) = -1.2 with 2% free electrons models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 2% interpolated. No constraint in ionization parameter.'
            print('Interpolation for the AGN a(ox) = -1.2 with 2% free electrons models is going to be used.')
            print('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
            res_NO = 0.125

    # AGN MODEL FOR alpha_OX = -1.2, efrac = 98% and logU in [-4.0, -0.5]
    elif sed == 3 and alpha == 2 and efrac == 2:
        file_lib = 'C17_AGN_alpha12_efrac98_CNfix.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file11:
            for line in file11:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 98%. No interpolation. No constraint in ionization parameter.'
            print('No interpolation for the AGN a(ox) = -1.2 with 98% free electrons models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 98% interpolated. No constraint in ionization parameter.'
            print('Interpolation for the AGN a(ox) = -1.2 with 98% free electrons models is going to be used.')
            print('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')

    # Different library
    elif sed == 4:
        file_lib = new_library
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file12:
            for line in file12:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'User file ' + file_lib + ' used as library for the models no interpolated'
            print('No interpolation for the library ' + file_lib)
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'User file ' + file_lib + ' used as library for the models interpolated'
            print('Interpolation for the library ' + file_lib)

    # Valuable columns of the files
    opt_lin = ['12logOH', 'logNO', 'logU', 'OII_3727', 'NeIII_3868', 'OIII_4363', 'OIII_5007', 'NII_6584', 'SII_671731']
    lin_opt_label = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'OII_3727', 'NeIII_3868', 'OIII_4363', 'OIII_5007',
                     'NII_6584', 'SII_6717,31']
    # POPSTAR MODEL
    if sed == 1:
        file_lib = 'C17_POPSTAR_1myr.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file6:
            for line in file6:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. No interpolation'
            print('No interpolation for the POPSTAR models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            print('')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. Interpolation'
            print('Interpolation for the POPSTAR models is going to be used.')
            print('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
            print('')
            res_NO = 0.125

    # BPASS MODEL
    elif sed == 2:
        file_lib = 'C17_BPASS_IMF135_mup300_1myr.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file7:
            for line in file7:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'BPASS a_IMF = 1.35, M_up = 300, age = 1Myr. No interpolation'
            print('No interpolation for the BPASS models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            print('')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'BPASS v.2.1, a_IMF = 1.35, M_up = 300, age = 1Myr. Interpolation'
            print('Interpolation for the BPASS  models is going to be used.')
            print('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
            print('')
            res_NO = 0.125


    # AGN MODEL FOR alpha_OX = -0.8, efrac = 2% and logU in [-4.0, -0.5]
    elif sed == 3 and alpha == 1 and efrac == 1:
        file_lib = 'C17_AGN_alpha08_efrac02_CNfix.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file8:
            for line in file8:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 2%. No interpolation. No constraint in ionization parameter.'
            print('No interpolation for the AGN a(ox) = -0.8 with 2% free electrons models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O. ')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 2% interpolated. No constraint in ionization parameter.'
            print('Interpolation for the AGN a(ox) = -0.8 with 2% free electrons models is going to be used.')
            print('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
            res_NO = 0.125

    # AGN MODEL FOR alpha_OX = -0.8, efrac = 98% and logU in [-4.0, -0.5]
    elif sed == 3 and alpha == 1 and efrac == 2:
        file_lib = 'C17_AGN_alpha08_efrac98_CNfix.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file9:
            for line in file9:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 98%. No interpolation. No constraint in ionization parameter.'
            print('No interpolation for the AGN a(ox) = -0.8 with 98% free electrons models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -0.8 and free electron fraction = 98% interpolated. No constraint in ionization parameter.'
            print('Interpolation for the AGN a(ox) = -0.8 with 98% free electrons models is going to be used.')
            print('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
            res_NO = 0.125

    # AGN MODEL FOR alpha_OX = -1.2, efrac = 2% and logU in [-4.0, -0.5]
    elif sed == 3 and alpha == 2 and efrac == 1:
        file_lib = 'C17_AGN_alpha12_efrac02_CNfix.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file10:
            for line in file10:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 2%. No interpolation. No constraint in ionization parameter.'
            print('No interpolation for the AGN a(ox) = -1.2 with 2% free electrons models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 2% interpolated. No constraint in ionization parameter.'
            print('Interpolation for the AGN a(ox) = -1.2 with 2% free electrons models is going to be used.')
            print('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
            res_NO = 0.125

    # AGN MODEL FOR alpha_OX = -1.2, efrac = 98% and logU in [-4.0, -0.5]
    elif sed == 3 and alpha == 2 and efrac == 2:
        file_lib = 'C17_AGN_alpha12_efrac98_CNfix.dat'
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file11:
            for line in file11:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 98%. No interpolation. No constraint in ionization parameter.'
            print('No interpolation for the AGN a(ox) = -1.2 with 98% free electrons models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -1.2 and free electron fraction = 98% interpolated. No constraint in ionization parameter.'
            print('Interpolation for the AGN a(ox) = -1.2 with 98% free electrons models is going to be used.')
            print('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')

    # Different library
    elif sed == 4:
        file_lib = HCm_folder
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib, 'r') as file12:
            for line in file12:
                if line[0] == '#':
                    n_comments += 1
        grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                                 skip_header=n_comments)
        if inter == 0:
            sed_type = 'User file ' + file_lib + ' used as library for the models no interpolated'
            print('No interpolation for the library ' + file_lib)
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'User file ' + file_lib + ' used as library for the models interpolated'
            print('Interpolation for the library ' + file_lib)

    # Valuable columns of the files
    opt_lin = ['12logOH', 'logNO', 'logU', 'OII_3727', 'NeIII_3868', 'OIII_4363', 'OIII_5007', 'NII_6584', 'SII_671731']
    lin_opt_label = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'OII_3727', 'NeIII_3868', 'OIII_4363', 'OIII_5007',
                     'NII_6584', 'SII_6717,31']

    ########################################
    ###### SORTING THE GRID OF MODELS ######
    ########################################

    print(' ')
    print('Sorting the grid of models')
    print(' ')

    index_OH_NO_U_sorted = []  # storing the correct order of the indexes

    # Sorting abundances 12+log(O/H)
    OH_values = grid_aux['12logOH']  # Oxygen abundances
    if len(OH_values) != 1:
        sorted_list_OH = sorted(range(len(OH_values)), key=OH_values.__getitem__)
    if len(OH_values) == 1:
        sorted_list_OH = [0]

    # Sorting abundance ratios log(N/O)
    OH_values_diff = list(set(OH_values[sorted_list_OH]))
    OH_values_diff.sort()  # It is necessary to sort again the list of different elements
    for OH_num in OH_values_diff:
        index_OH_fix = np.where(OH_values == OH_num)[0]  # Index(es) for a particular abundance 12+log(O/H)
        NO_values = grid_aux['logNO'][index_OH_fix]
        if len(NO_values) != 1:
            sorted_list_NO = sorted(range(len(NO_values)), key=NO_values.__getitem__)
        if len(NO_values) == 1:
            sorted_list_NO = [0]
        NO_values_diff = list(set(NO_values[sorted_list_NO]))
        NO_values_diff.sort()  # It s necessary to sort again the list of different elements
        for NO_num in NO_values_diff:
            index_OH_NO_fix = np.where(NO_values == NO_num)[
                0]  # Index(es) for particular abundances 12+log(O/H) and log(N/O)
            # Sorting ionization parameters
            U_values = grid_aux['logU'][index_OH_fix[index_OH_NO_fix]]
            if len(U_values) != 1:
                sorted_list_U = sorted(range(len(U_values)), key=U_values.__getitem__)
            if len(U_values) == 1:
                sorted_list_U = [0]
            index_OH_NO_U = index_OH_fix[index_OH_NO_fix[sorted_list_U]]  # Sorted index(es) for U at fixed O/H and N/O
            for index_sort in index_OH_NO_U:
                index_OH_NO_U_sorted.append(index_sort)  # Adding index in the correct order

    # Generating new library file
    list_comments = []  # Storing comments in the file:
    with open('Libraries_opt/' + file_lib, 'r') as file_aux:
        for line in file_aux:
            if line[0] == '#':
                list_comments.append(line)

    # Storing columns:
    lin_opt_col = []
    # Retrieving each column of the grid
    for label in opt_lin:
        aux_col = grid_aux[label].tolist()
        lin_opt_col.append(aux_col)

    # Comments
    grid_to_write = open('Libraries_opt/' + file_lib, 'w')
    for line_com in list_comments:
        grid_to_write.write(line_com)
    # Header line
    label_line = '{:15} '.format(lin_opt_label[0].replace(' ', ''))
    for ind in range(1, len(lin_opt_label) - 1):
        label_line += '\t {:15} '.format(lin_opt_label[ind].replace(' ', ''))
    label_line += '\t {:15}\n'.format(lin_opt_label[-1].replace(' ', ''))
    grid_to_write.write(label_line)
    # Values:
    for ind_val in index_OH_NO_U_sorted:
        val_line = '{:7.7f} '.format(lin_opt_col[0][ind_val])
        for ind2 in range(1, len(lin_opt_label) - 1):
            val_line += '\t {:7.7f} '.format(lin_opt_col[ind2][ind_val])
        val_line += '\t {:7.7f}\n'.format(lin_opt_col[-1][ind_val])
        grid_to_write.write(val_line)
    grid_to_write.close()

    # Opening sorted grid of models
    n_comments = 0
    with open('Libraries_opt/' + file_lib, 'r') as file12:
        for line in file12:
            if line[0] == '#':
                n_comments += 1
    grid_aux = np.genfromtxt('Libraries_opt/' + file_lib, dtype=None, names=True, encoding='ascii',
                             skip_header=n_comments)

    ################################################
    ###### CONSTRAINTS FOR THE GRID OF MODELS ######
    ################################################

    # Reading constraints and creating library with constraints
    print('')
    print(
        'Select a file with the constraints to be used to limit the grid of models when the measurement of a quantity is impossible without any relation.')
    print('')

    # Displayig options for the user
    print('')
    question = True
    while question:
        print('-------------------------------------------------')
        print('Default constraints')
        print('-------------------')
        print('(1) Constraints for Star-Forming Galaxies')
        print('(2) Constraints for Extreme Emission Line Galaxies')
        print('(3) Constraints for AGNs (no restriction in the ionization parameter)')
        print('(4) Constraints for high ionization AGNs (log(U) > -2.5)')
        print('(5) Constraints for low ionization AGNs (log(U) < -2.5)')
        print('')
        print('Other constraints')
        print('-----------------')
        print('(6) Different constraint file')
        print('-------------------------------------------------')
        # if int(sys.version[0]) < 3:
        #     const = raw_input('Choose constraint for the grids: ')
        # else:
        #     const = input('Choose constraint for the grids: ')
        # if const == '1' or const == '2' or const == '3' or const == '4' or const == '5' or const == '6': question = False

    # Particular file introduced by the user
    if const == '6':
        question = True
        while question:
            print(
                'Introduce name of the file containing the constraints for the grids. It must be located in the folder "Constraints".')
            print('')
            if int(sys.version[0]) < 3:
                new_const = input('Name of file: ')
            else:
                new_const = input('Name of file: ')

            # Searching for the file
            try:
                # Counting comments:
                n_comments = 0
                with open('Constraints/' + new_const, 'r') as file13:
                    for line in file13:
                        if line[0] == '#':
                            n_comments += 1
                const_user = np.genfromtxt('Constraints/' + new_const, dtype=None, names=True, encoding='ascii',
                                           skip_header=n_comments)
                print('')
                print('Loading constraint file ' + new_const + '. Checking correct format of the file.')
                question = False
            except:
                print('')
                print('File was not found in folder "Constraints" or file does not exist.')
        question = True
        while question:
            try:
                # Counting comments:
                n_comments = 0
                with open('Constraints/' + new_const, 'r') as file14:
                    for line in file14:
                        if line[0] == '#':
                            n_comments += 1
                const_user = np.genfromtxt('Constraints/' + new_const, dtype=None, names=True, encoding='ascii',
                                           skip_header=n_comments)
                # Checking correct format:
                # Counting comments:
                n_comments = 0
                with open('Constraints/template_OH.dat', 'r') as file15:
                    for line in file15:
                        if line[0] == '#':
                            n_comments += 1
                auxiliar_labels = np.genfromtxt('Constraints/template_OH.dat', dtype=None, names=True, encoding='ascii',
                                                skip_header=n_comments).dtype.names
                missing_labels = []
                for label in auxiliar_labels:
                    if label in const_user.dtype.names:
                        continue
                    else:
                        missing_labels.append(label)
                # Displaying message for the user:
                print('Succesfully reading of the file')
                if len(missing_labels) == 0:
                    print('File presents the correct format')
                    question = False
                else:
                    print('File does not present the correct format. The following columns are missing:')
                    for need_label in missing_labels:
                        print('- ' + need_label)
                    print('More details on the correct format for the library are found in readme file.')
                    print('')
                    print('Reintroduce name of the file with fixed format.')
                    print('')
                    if int(sys.version[0]) < 3:
                        new_const = input('Name of file: ')
                    else:
                        new_const = input('Name of file: ')
            except:
                print('Something went wrong while reading file. Please, reintroduce name of the file.')
                print('')
                if int(sys.version[0]) < 3:
                    new_const = input('Name of file: ')
                else:
                    new_const = input('Name of file: ')

    # Generation of grids with constrained laws:
    if const == '1' or const == '2' or const == '3' or const == '6':
        # First grid does not change
        grid1 = grid_aux
        file_lib_2 = file_lib
    elif const == '4' or const == '5':
        lin_opt_agn = []
        # The initial grid need to be constrained in the ionization parameter
        if const == '4':
            U_max = 0.0
            U_min = -2.5
            tag = 'high'
        if const == '5':
            U_max = -2.5
            U_min = -4.0
            tag = 'low'
        # Retrieving each column of the grid
        for label in opt_lin:
            aux_col = grid_aux[label].tolist()
            lin_opt_agn.append(aux_col)
        # Creation of the grid
        file_lib_2 = '.'.join(file_lib.split('.')[0:-1]) + '_' + tag + '.' + file_lib.split('.')[-1]
        file_open = open('Libraries_opt/' + file_lib_2, 'w')
        file_open.write('#Library constrained for ' + tag + ' ionization AGNs\n')
        # Header line
        label_line = '{:15} '.format(lin_opt_label[0].replace(' ', ''))
        for ind in range(1, len(lin_opt_label) - 1):
            label_line += '\t {:15} '.format(lin_opt_label[ind].replace(' ', ''))
        label_line += '\t {:15}\n'.format(lin_opt_label[-1].replace(' ', ''))
        file_open.write(label_line)
        # Values:
        for ind_val in range(0, len(lin_opt_agn[0])):
            if lin_opt_agn[2][ind_val] <= U_max and lin_opt_agn[2][ind_val] >= U_min:
                val_line = '{:7.7f} '.format(lin_opt_agn[0][ind_val])
                for ind2 in range(1, len(lin_opt_label) - 1):
                    val_line += '\t {:7.7f} '.format(lin_opt_agn[ind2][ind_val])
                val_line += '\t {:7.7f}\n'.format(lin_opt_agn[-1][ind_val])
                file_open.write(val_line)
        file_open.close()
        # Counting comments:
        n_comments = 0
        with open('Libraries_opt/' + file_lib_2, 'r') as file:
            for line in file:
                if line[0] == '#':
                    n_comments += 1
        grid1 = np.genfromtxt('Libraries_opt/' + file_lib_2, dtype=None, names=True, encoding='ascii',
                              skip_header=n_comments)

    # Generating libraries for the constraints in the files
    if const == '1':  # Star-Forming Galaxies
        name_const = 'Constraints/template_OH.dat'
        const_file = 'template_OH.dat'
        n_comments = 0
        with open(name_const, 'r') as file16:
            for line in file16:
                if line[0] == '#':
                    n_comments += 1
        const_data = np.genfromtxt(name_const, dtype=None, names=True, encoding='ascii', skip_header=n_comments)
    if const == '2':
        name_const = 'Constraints/template_OH_eelg.dat'
        const_file = 'template_OH_eelg.dat'
        n_comments = 0
        with open(name_const, 'r') as file17:
            for line in file17:
                if line[0] == '#':
                    n_comments += 1
        const_data = np.genfromtxt(name_const, dtype=None, names=True, encoding='ascii', skip_header=n_comments)
    if const == '3':
        name_const = 'Constraints/template_OH_agn.dat'
        const_file = 'template_OH_agn.dat'
        n_comments = 0
        with open(name_const, 'r') as file18:
            for line in file18:
                if line[0] == '#':
                    n_comments += 1
        const_data = np.genfromtxt(name_const, dtype=None, names=True, encoding='ascii', skip_header=n_comments)
    if const == '4':
        name_const = 'Constraints/template_OH_agn_high.dat'
        const_file = 'template_OH_agn_high.dat'
        n_comments = 0
        with open(name_const, 'r') as file19:
            for line in file19:
                if line[0] == '#':
                    n_comments += 1
        const_data = np.genfromtxt(name_const, dtype=None, names=True, encoding='ascii', skip_header=n_comments)
    if const == '5':
        name_const = 'Constraints/template_OH_agn_low.dat'
        const_file = 'template_OH_agn_low.dat'
        n_comments = 0
        with open(name_const, 'r') as file20:
            for line in file20:
                if line[0] == '#':
                    n_comments += 1
        const_data = np.genfromtxt(name_const, dtype=None, names=True, encoding='ascii', skip_header=n_comments)
    if const == '6':
        const_file = new_const
        name_const = 'Constraints/' + new_const
        n_comments = 0
        with open(name_const, 'r') as file21:
            for line in file21:
                if line[0] == '#':
                    n_comments += 1
        const_data = np.genfromtxt(name_const, dtype=None, names=True, encoding='ascii', skip_header=n_comments)

    print('')

    return output

    return

def epm_HII_CHI_mistry_orig(input00, output_file, n, sed, inter, HCm_folder=None):

    if HCm_folder is None:
        HCm_folder = Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v5.1/')

    print(' ---------------------------------------------------------------------')
    print('This is HII-CHI-mistry v. 5.1')
    print(' See Perez-Montero, E. (2014) for details')
    print(' Insert the name of your input text file with all or some of the following columns:')
    print(' 3727 [OII], 3868 [NeIII], 4363 [OIII], 5007 [OIII], 6584 [NII], 6725 [SII]')
    print('with their corresponding labels and errors in adjacent columns')
    print('---------------------------------------------------------------------')

    # # Input file reading
    # try:
    #     input0 = np.genfromtxt(input00, dtype=None, names=True, encoding='ascii')
    #     print('The input file is:' + input00)
    # except:
    #     print('Input file error: It does not exist or has wrong format')
    #
    #     sys.exit
    #
    # if input0.size == 1:
    #     input1 = np.stack((input0, input0))
    # else:
    #     input1 = input0
    input0 = input00
    input1 = input00

    print('-------------------------------------------------')
    print('(1) POPSTAR with Chabrier IMF, age = 1 Myr')
    print('(2) BPASS v.2.1 a_IMF = 1.35, Mup = 300, age = 1Myr')
    print('(3) AGN, double component, a(OX) = -0.8, a(UV) = -1.0')
    print('(4) AGN, double component, a(OX) = -1.2, a(UV) = -1.0')
    print(f'-- SED selected: {sed}')

    print('-------------------------------------------------')
    print('Choose models [0] No interpolated [1] Interpolated:')
    print(f'-- Mode selected {inter}')

    sed = int(sed)
    inter = int(inter)

    grid_folder = Path(HCm_folder)

    if sed == 1:
        grid1 = np.loadtxt(grid_folder/'C17_popstar_v5.0.dat')
        grid2 = np.loadtxt(grid_folder/'C17_popstar_logU_adapted_emp_v5.0.dat')
        grid3 = np.loadtxt(grid_folder/'C17_popstar_logU-NO_adapted_emp_v5.0.dat')
        if inter == 0:
            sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. No interpolation'
            print('No interpolation for the POPSTAR models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            print('')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'POPSTAR, age = 1 Myr, Chabrier IMF. interpolation'
            print('Interpolation for the POPSTAR models is going to be used.')
            print('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
            print('')
            res_NO = 0.0125
    elif sed == 2:
        grid1 = np.loadtxt(grid_folder/'C17_bpass_v5.0.dat')
        grid2 = np.loadtxt(grid_folder/'C17_bpass_logU_adapted_emp_v5.0.dat')
        grid3 = np.loadtxt(grid_folder/'C17_bpass_logU-NO_adapted_emp_v5.0.dat')
        if inter == 0:
            sed_type = 'BPASS a_IMF = 1.35, M_up = 300, age = 1Myr. No interpolation'
            print('No interpolation for theBPASS models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            print('')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'BPASS v.2.1, a_IMF = 1.35, M_up = 300, age = 1Myr. Interpolation'
            print('Interpolation for theBPASS  models is going to be used.')
            print('The grid has a resolution of 0.01dex for O/H and 0.0125dex for N/O')
            print('')
            res_NO = 0.0125
    elif sed == 3:
        grid1 = np.loadtxt(grid_folder/'C17_agn_a08_v5.0.dat')
        grid2 = np.loadtxt(grid_folder/'C17_agn_a08_v5.0.dat')
        grid3 = np.loadtxt(grid_folder/'C17_agn_a08_NO_adapted_emp_v5.0.dat')
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -0.8. No interpolation'
            print('No interpolation for the AGN a(ox) = -0.8 models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            print('')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -0.8. Interpolation'
            print('Interpolation for the AGN a(ox) = -0.8 models is going to be used.')
            print('The grid has a resolution of 0.01dex for O/H and 0.0125 dex for N/O')
            print('')
            res_NO = 0.0125
    elif sed == 4:
        grid1 = np.loadtxt(grid_folder/'C17_agn_a12_v5.0.dat')
        grid2 = np.loadtxt(grid_folder/'C17_agn_a12_v5.0.dat')
        grid3 = np.loadtxt(grid_folder/'C17_agn_a12_NO_adapted_emp_v5.0.dat')
        if inter == 0:
            sed_type = 'Double composite AGN, a(OX) = -1.2. No interpolation'
            print('No interpolation for the AGN a(ox) = -1.2 models is going to be used.')
            print('The grid has a resolution of 0.1dex for O/H and 0.125dex for N/O')
            print('')
            res_NO = 0.125
        elif inter == 1:
            sed_type = 'Double composite AGN, a(OX) = -1.2. Interpolation'
            print('Interpolation for the AGN a(ox) = -1.2 models is going to be used.')
            print('The grid has a resolution of 0.01 dex for O/H and 0.0125 dex for N/O')
            print('')
            res_NO = 0.125

    grids = []
    OHffs = []
    eOHffs = []
    NOffs = []
    eNOffs = []
    logUffs = []
    elogUffs = []

    Label_ID = False
    Label_OII = False
    Label_eOII = False
    Label_NeIII = False
    Label_eNeIII = False
    Label_OIIIa = False
    Label_eOIIIa = False
    Label_OIII_4959 = False
    Label_eOIII_4959 = False
    Label_OIII_5007 = False
    Label_eOIII_5007 = False
    Label_NII = False
    Label_eNII = False
    Label_SII = False
    Label_eSII = False
    Label_SII_6716 = False
    Label_eSII_6716 = False
    Label_SII_6731 = False
    Label_eSII_6731 = False

    for col in range(0, len(input1.dtype.names), 1):
        if input1.dtype.names[col] == 'ID':
            Label_ID = True
        if input1.dtype.names[col] == 'OII_3727':
            Label_OII = True
        if input1.dtype.names[col] == 'eOII_3727':
            Label_eOII = True
        if input1.dtype.names[col] == 'NeIII_3868':
            Label_NeIII = True
        if input1.dtype.names[col] == 'eNeIII_3868':
            Label_eNeIII = True
        if input1.dtype.names[col] == 'OIII_4363':
            Label_OIIIa = True
        if input1.dtype.names[col] == 'eOIIIa_4363':
            Label_eOIIIa = True
        if input1.dtype.names[col] == 'OIII_4959':
            Label_OIII_4959 = True
        if input1.dtype.names[col] == 'eOIII_4959':
            Label_eOIII_4959 = True
        if input1.dtype.names[col] == 'OIII_5007':
            Label_OIII_5007 = True
        if input1.dtype.names[col] == 'eOIII_5007':
            Label_eOIII_5007 = True
        if input1.dtype.names[col] == 'NII_6584':
            Label_NII = True
        if input1.dtype.names[col] == 'eNII_6584':
            Label_eNII = True
        if input1.dtype.names[col] == 'SII_6725':
            Label_SII = True
        if input1.dtype.names[col] == 'eSII_6725':
            Label_eSII = True
        if input1.dtype.names[col] == 'SII_6716':
            Label_SII_6716 = True
        if input1.dtype.names[col] == 'SII_6731':
            Label_SII_6731 = True
        if input1.dtype.names[col] == 'eSII_6716':
            Label_eSII_6716 = True
        if input1.dtype.names[col] == 'eSII_6731':
            Label_eSII_6731 = True

    if Label_ID == False:
        Names = np.arange(1, input1.size + 1, 1)
    else:
        Names = input1['ID']
    if Label_OII == False:
        OII_3727 = np.zeros(input1.size)
    else:
        OII_3727 = input1['OII_3727']
    if Label_eOII == False:
        eOII_3727 = np.zeros(input1.size)
    else:
        eOII_3727 = input1['eOII_3727']
    if Label_NeIII == False:
        NeIII_3868 = np.zeros(input1.size)
    else:
        NeIII_3868 = input1['NeIII_3868']
    if Label_eNeIII == False:
        eNeIII_3868 = np.zeros(input1.size)
    else:
        eNeIII_3868 = input1['eNeIII_3868']
    if Label_OIIIa == False:
        OIII_4363 = np.zeros(input1.size)
    else:
        OIII_4363 = input1['OIII_4363']
    if Label_eOIIIa == False:
        eOIII_4363 = np.zeros(input1.size)
    else:
        eOIII_4363 = input1['eOIII_4363']
    if Label_OIII_4959 == False and Label_OIII_5007 == False:
        OIII_5007 = np.zeros(input1.size)
    elif Label_OIII_4959 == False and Label_OIII_5007 == True:
        OIII_5007 = input1['OIII_5007']
    elif Label_OIII_4959 == True and Label_OIII_5007 == False:
        OIII_5007 = 3 * input1['eOIII_4959']
    else:
        OIII_5007 = (input1['OIII_5007'] + input1['OIII_4959']) / 1.33
    if Label_eOIII_4959 == False and Label_eOIII_5007 == False:
        eOIII_5007 = np.zeros(input1.size)
    elif Label_eOIII_4959 == False and Label_eOIII_5007 == True:
        eOIII_5007 = input1['eOIII_5007']
    elif Label_eOIII_4959 == True and Label_eOIII_5007 == False:
        eOIII_5007 = 4 * input1['eOIII_4959']
    else:
        eOIII_5007 = (input1['eOIII_5007'] + input1['eOIII_4959']) / 1.33
    if Label_NII == False:
        NII_6584 = np.zeros(input1.size)
    else:
        NII_6584 = input1['NII_6584']
    if Label_eNII == False:
        eNII_6584 = np.zeros(input1.size)
    else:
        eNII_6584 = input1['eNII_6584']
    if Label_SII == False and (Label_SII_6716 == False or Label_SII_6731 == False):
        SII_6725 = np.zeros(input1.size)
    elif Label_SII == True:
        SII_6725 = input1['SII_6725']
    else:
        SII_6725 = input1['SII_6716'] + input1['SII_6731']
    if Label_eSII == False and (Label_eSII_6716 == False or Label_eSII_6731 == False):
        eSII_6725 = np.zeros(input1.size)
    elif Label_eSII == True:
        eSII_6725 = input1['eSII_6725']
    else:
        eSII_6725 = input1['eSII_6716'] + input1['eSII_6731']

    output = np.zeros(input1.size,
                      dtype=[('ID', 'U12'), ('OII_3727', float), ('eOII_3727', float), ('NeIII_3868', float),
                             ('eNeIII_3868', float), ('OIII_4363', float), ('eOIII_4363', float), ('OIII_5007', float),
                             ('eOIII_5007', float), ('NII_6584', float), ('eNII_6584', float), ('SII_6725', float),
                             ('eSII_6725', float), ('grid', int), ('OH', float), ('eOH', float), ('NO', float),
                             ('eNO', float), ('logU', float), ('elogU', float)])

    output['ID'] = Names
    output['OII_3727'] = OII_3727
    output['eOII_3727'] = eOII_3727
    output['NeIII_3868'] = NeIII_3868
    output['eNeIII_3868'] = eNeIII_3868
    output['OIII_4363'] = OIII_4363
    output['eOIII_4363'] = eOIII_4363
    output['OIII_5007'] = OIII_5007
    output['eOIII_5007'] = eOIII_5007
    output['NII_6584'] = NII_6584
    output['eNII_6584'] = eNII_6584
    output['SII_6725'] = SII_6725
    output['eSII_6725'] = eSII_6725



    print('Reading grids ....')
    print('')
    print('')
    print('----------------------------------------------------------------')
    print('(%)   ID    Grid  12+log(O/H)  log(N/O)    log(U)')
    print('-----------------------------------------------------------------')

    count = 0
    for tab in range(0, input1.size, 1):
        count = count + 1

        OH_mc = []
        NO_mc = []
        logU_mc = []
        OHe_mc = []
        NOe_mc = []
        logUe_mc = []

        # Selection of grid

        if OIII_4363[tab] > 0 and OIII_5007[tab] > 0:
            grid = grid1
            grid_type = 1
            grids.append(1)
        elif NII_6584[tab] > 0 and (OII_3727[tab] > 0 or SII_6725[tab] > 0):
            grid = grid2
            grid_type = 2
            grids.append(2)
        else:
            grid = grid3
            grid_type = 3
            grids.append(3)

        # Calculation of N/O

        if NII_6584[tab] == 0 or (OII_3727[tab] == 0 and SII_6725[tab] == 0):
            NOff = -10
            eNOff = 0
        else:
            for monte in range(0, n, 1):
                NO_p = 0
                den_NO = 0
                NO_e = 0
                den_NO_e = 0
                tol_max = 1e2

                OII_3727_obs = 0
                if OII_3727[tab] == 0:
                    OII_3727_obs = 0
                else:
                    while OII_3727_obs <= 0:
                        OII_3727_obs = np.random.normal(OII_3727[tab], eOII_3727[tab] + 1e-3)
                OIII_4363_obs = 0
                if OIII_4363[tab] == 0:
                    OIII_4363_obs = 0
                else:
                    while OIII_4363_obs <= 0:
                        OIII_4363_obs = np.random.normal(OIII_4363[tab], eOIII_4363[tab] + 1e-3)
                OIII_5007_obs = 0
                if OIII_5007[tab] == 0:
                    OIII_5007_obs = 0
                else:
                    while OIII_5007_obs <= 0:
                        OIII_5007_obs = np.random.normal(OIII_5007[tab], eOIII_5007[tab] + 1e-3)
                if OIII_4363_obs == 0 or OIII_5007_obs == 0:
                    ROIII_obs = 0
                else:
                    ROIII_obs = OIII_5007_obs / OIII_4363_obs
                NII_6584_obs = 0
                if NII_6584[tab] == 0:
                    NII_6584_obs = 0
                else:
                    while NII_6584_obs <= 0:
                        NII_6584_obs = np.random.normal(NII_6584[tab], eNII_6584[tab] + 1e-3)
                SII_6725_obs = 0
                if SII_6725[tab] == 0:
                    SII_6725_obs = 0
                else:
                    while SII_6725_obs <= 0:
                        SII_6725_obs = np.random.normal(SII_6725[tab], eSII_6725[tab] + 1e-3)
                    if SII_6725_obs <= 0: SII_6725_obs = 0
                if NII_6584_obs == 0 or OII_3727_obs == 0:
                    N2O2_obs = -10
                else:
                    N2O2_obs = np.log10(NII_6584_obs / OII_3727_obs)
                if NII_6584_obs == 0 or SII_6725_obs == 0:
                    N2S2_obs = -10
                else:
                    N2S2_obs = np.log10(NII_6584_obs / SII_6725_obs)

                CHI_ROIII = 0
                CHI_N2O2 = 0
                CHI_N2S2 = 0
                CHI_NO = 0

                for index in grid:
                    if ROIII_obs == 0:
                        CHI_ROIII = 0
                    elif index[5] == 0:
                        CHI_ROIII = tol_max
                    else:
                        CHI_ROIII = (index[6] / index[5] - ROIII_obs) ** 2 / (index[6] / index[5])
                    if N2O2_obs == -10:
                        CHI_N2O2 = 0
                    elif index[3] == 0 or index[7] == 0:
                        CHI_N2O2 = tol_max
                    else:
                        CHI_N2O2 = (np.log10(index[7] / index[3]) - N2O2_obs) ** 2 / (
                            abs(np.log10(index[7] / index[3]) + 1e-3))
                    if N2S2_obs == -10:
                        CHI_N2S2 = 0
                    elif index[7] == 0 or index[8] == 0:
                        CHI_N2S2 = tol_max
                    else:
                        CHI_N2S2 = (np.log10(index[7] / index[8]) - N2S2_obs) ** 2 / (
                            abs(np.log10(index[7] / index[8]) + 1e-3))
                    CHI_NO = (CHI_ROIII ** 2 + CHI_N2O2 ** 2 + CHI_N2S2 ** 2) ** 0.5

                    if CHI_NO == 0:
                        NO_p = NO_p
                        den_NO = den_NO
                    else:
                        NO_p = index[1] / (CHI_NO) + NO_p
                        den_NO = 1 / (CHI_NO) + den_NO

                NO = NO_p / den_NO

                # Calculation of N/O error

                CHI_ROIII = 0
                CHI_N2O2 = 0
                CHI_N2S2 = 0
                CHI_NO = 0

                for index in grid:
                    if ROIII_obs == 0:
                        CHI_ROIII = 0
                    elif index[5] == 0:
                        CHI_ROIII = tol_max
                    else:
                        CHI_ROIII = (index[6] / index[5] - ROIII_obs) ** 2 / (index[6] / index[5])
                    if N2O2_obs == -10:
                        CHI_N2O2 = 0
                    elif index[3] == 0 or index[7] == 0:
                        CHI_N2O2 = tol_max
                    else:
                        CHI_N2O2 = (np.log10(index[7] / index[3]) - N2O2_obs) ** 2 / (
                            abs(np.log10(index[7] / index[3]) + 1e-3))
                    if N2S2_obs == -10:
                        CHI_N2S2 = 0
                    elif index[7] == 0 or index[8] == 0:
                        CHI_N2S2 = tol_max
                    else:
                        CHI_N2S2 = (np.log10(index[7] / index[8]) - N2S2_obs) ** 2 / (
                            abs(np.log10(index[7] / index[8]) + 1e-3))

                    CHI_NO = (CHI_ROIII ** 2 + CHI_N2O2 ** 2 + CHI_N2S2 ** 2) ** 0.5

                    if CHI_NO == 0:
                        NO_e = NO_e
                        den_NO_e = den_NO_e
                    else:
                        NO_e = (index[1] - NO) ** 2 / (CHI_NO) + NO_e
                        den_NO_e = 1 / (CHI_NO) + den_NO_e

                eNO = NO_e / den_NO_e

                # Iterations for the interpolation mode

                if inter == 0 or NO == -10:
                    NOf = NO
                elif inter == 1:
                    igrid = grid[np.lexsort((grid[:, 0], grid[:, 2]))]
                    igrid = interpolate_orig(igrid, 1, NO - eNO - 0.125, NO + eNO + 0.125, 10)

                    CHI_ROIII = 0
                    CHI_N2O2 = 0
                    CHI_N2S2 = 0
                    CHI_NO = 0
                    NO_p = 0
                    den_NO = 0

                    for index in igrid:
                        if ROIII_obs == 0:
                            CHI_ROIII = 0
                        elif index[5] == 0:
                            CHI_ROIII = tol_max
                        else:
                            CHI_ROIII = (index[6] / index[5] - ROIII_obs) ** 2 / (index[6] / index[5])
                            if OIII_5007_obs == 0:
                                CHI_OIII = 0
                            elif index[6] == 0:
                                CHI_OIII = tol_max
                            else:
                                CHI_OIII = (index[6] - OIII_5007_obs) ** 2 / index[6]
                            if OII_3727_obs == 0:
                                CHI_OII = 0
                            elif index[3] == 0:
                                CHI_OII = tol_max
                            else:
                                CHI_OII = (index[3] - OII_3727_obs) ** 2 / index[3]
                        if N2O2_obs == -10:
                            CHI_N2O2 = 0
                        elif index[3] == 0 or index[7] == 0:
                            CHI_N2O2 = tol_max
                        else:
                            CHI_N2O2 = (np.log10(index[7] / index[3]) - N2O2_obs) ** 2 / (
                                abs(np.log10(index[7] / index[3]) + 1e-3))
                        if N2S2_obs == -10:
                            CHI_N2S2 = 0
                        elif index[7] == 0 or index[8] == 0:
                            CHI_N2S2 = tol_max
                        else:
                            CHI_N2S2 = (np.log10(index[7] / index[8]) - N2S2_obs) ** 2 / (
                                abs(np.log10(index[7] / index[8]) + 1e-3))

                        CHI_NO = (CHI_ROIII ** 2 + CHI_N2O2 ** 2 + CHI_N2S2 ** 2) ** 0.5
                        if CHI_NO == 0:
                            NO_p = NO_p
                            den_NO = den_NO
                        else:
                            NO_p = index[1] / CHI_NO + NO_p
                            den_NO = 1 / CHI_NO + den_NO

                    NOf = NO_p / den_NO

                NO_mc.append(NOf)
                NOe_mc.append(eNO)

            NOff = np.mean(NO_mc)
            eNOff = (np.std(NO_mc) ** 2 + np.mean(NOe_mc) ** 2) ** 0.5

        # Creation of a constrained grid on N/O
        if NOff == -10:
            grid_c = grid
        else:
            grid_mac = []
            for index in grid:
                if np.abs(index[1] - NOff) > np.abs(eNOff + res_NO):
                    continue
                else:
                    grid_mac.append(index[0])
                    grid_mac.append(index[1])
                    grid_mac.append(index[2])
                    grid_mac.append(index[3])
                    grid_mac.append(index[4])
                    grid_mac.append(index[5])
                    grid_mac.append(index[6])
                    grid_mac.append(index[7])
                    grid_mac.append(index[8])
                grid_c = np.reshape(grid_mac, (int(len(grid_mac) / 9), 9))

        # Calculation of O/H and logU
        for monte in range(0, n, 1):

            OH_p = 0
            logU_p = 0
            den_OH = 0
            OH_e = 0
            logU_e = 0
            den_OH_e = 0
            tol_max = 1e2

            OII_3727_obs = 0
            if OII_3727[tab] == 0:
                OII_3727_obs = 0
            else:
                while OII_3727_obs <= 0:
                    OII_3727_obs = np.random.normal(OII_3727[tab], eOII_3727[tab] + 1e-3)
            NeIII_3868_obs = 0
            if NeIII_3868[tab] == 0:
                NeIII_3868_obs = 0
            else:
                while NeIII_3868_obs <= 0:
                    NeIII_3868_obs = np.random.normal(NeIII_3868[tab], eNeIII_3868[tab] + 1e-3)
            OIII_4363_obs = 0
            if OIII_4363[tab] == 0:
                OIII_4363_obs = 0
            else:
                while OIII_4363_obs <= 0:
                    OIII_4363_obs = np.random.normal(OIII_4363[tab], eOIII_4363[tab] + 1e-3)
            OIII_5007_obs = 0
            if OIII_5007[tab] == 0:
                OIII_5007_obs = 0
            else:
                while OIII_5007_obs <= 0:
                    OIII_5007_obs = np.random.normal(OIII_5007[tab], eOIII_5007[tab] + 1e-3)
            if OIII_4363_obs == 0 or OIII_5007_obs == 0:
                ROIII_obs = 0
            else:
                ROIII_obs = OIII_5007_obs / OIII_4363_obs
            NII_6584_obs = 0
            if NII_6584[tab] == 0:
                NII_6584_obs = 0
            else:
                while NII_6584_obs <= 0:
                    NII_6584_obs = np.random.normal(NII_6584[tab], eNII_6584[tab] + 1e-3)
            SII_6725_obs = 0
            if SII_6725[tab] == 0:
                SII_6725_obs = 0
            else:
                while SII_6725_obs <= 0:
                    SII_6725_obs = np.random.normal(SII_6725[tab], eSII_6725[tab] + 1e-3)
            if OII_3727_obs == 0 or OIII_5007_obs == 0:
                O2O3_obs = 0
                R23_obs = -10
            else:
                R23_obs = np.log10(OII_3727_obs + OIII_5007_obs)
                O2O3_obs = (OII_3727_obs / OIII_5007_obs)
            if OII_3727_obs == 0 or NeIII_3868_obs == 0:
                O2Ne3_obs = 0
                R2Ne3_obs = -10
            else:
                O2Ne3_obs = (OII_3727_obs / NeIII_3868_obs)
                R2Ne3_obs = np.log10(OII_3727_obs + NeIII_3868_obs)
            if OIII_5007_obs == 0 or NII_6584_obs == 0:
                O3N2_obs = -10
            else:
                O3N2_obs = np.log10(OIII_5007_obs / NII_6584_obs)
            if OIII_5007_obs == 0 or SII_6725_obs == 0:
                O3S2_obs = -10
            else:
                O3S2_obs = np.log10(OIII_5007_obs / SII_6725_obs)

            if R23_obs == -10 and NII_6584_obs == 0 and ROIII_obs == 0 and R2Ne3_obs == -10 and O3S2_obs == -10:
                OH = 0
                logU = 0
            else:
                CHI_ROIII = 0
                CHI_NII = 0
                CHI_OIII = 0
                CHI_OII = 0
                CHI_O2O3 = 0
                CHI_R23 = 0
                CHI_O2Ne3 = 0
                CHI_R2Ne3 = 0
                CHI_O3N2 = 0
                CHI_O3S2 = 0
                CHI_OH = 0
                for index in grid_c:
                    if ROIII_obs == 0:
                        CHI_ROIII = 0
                    elif index[5] == 0:
                        CHI_ROIII = tol_max
                    else:
                        CHI_ROIII = (index[6] / index[5] - ROIII_obs) ** 2 / (index[6] / index[5])
                    if OIII_5007_obs == 0:
                        CHI_OIII = 0
                    elif index[6] == 0:
                        CHI_OIII = tol_max
                    else:
                        CHI_OIII = (index[6] - OIII_5007_obs) ** 2 / index[6]
                    if OII_3727_obs == 0:
                        CHI_OII = 0
                    elif index[3] == 0:
                        CHI_OII = tol_max
                    else:
                        CHI_OII = (index[3] - OII_3727_obs) ** 2 / index[3]
                    if NII_6584_obs == 0:
                        CHI_NII = 0
                    elif index[7] == 0:
                        CHI_NII = tol_max
                    else:
                        CHI_NII = (index[7] - NII_6584_obs) ** 2 / index[7]
                    if OII_3727_obs == 0 or OIII_5007_obs == 0:
                        CHI_O2O3 = 0
                        CHI_R23 = 0
                    elif index[3] == 0 or index[6] == 0:
                        CHI_O2O3 = tol_max
                        CHI_R23 = tol_max
                    else:
                        CHI_O2O3 = (index[3] / index[6] - O2O3_obs) ** 2 / (index[3] / index[6])
                        CHI_R23 = (np.log10(index[3] + index[6]) - R23_obs) ** 2 / (
                            np.abs(np.log10(index[3] + index[6] + 1e-3)))
                    if OII_3727_obs == 0 or NeIII_3868_obs == 0:
                        CHI_O2Ne3 = 0
                        CHI_R2Ne3 = 0
                    elif index[3] == 0 or index[4] == 0:
                        CHI_O2Ne3 = tol_max
                        CHI_R2Ne3 = tol_max
                    else:
                        CHI_O2Ne3 = (index[3] / index[4] - O2Ne3_obs) ** 2 / (index[3] / index[4])
                        CHI_R2Ne3 = (np.log10(index[3] + index[4]) - R2Ne3_obs) ** 2 / (
                            np.abs(np.log10(index[3] + index[4] + 1e-3)))
                    if OIII_5007_obs == 0 or NII_6584_obs == 0:
                        CHI_O3N2 = 0
                    elif index[6] == 0 or index[7] == 0:
                        CHI_O3N2 = tol_max
                    else:
                        CHI_O3N2 = (np.log10(index[6] / index[7]) - O3N2_obs) ** 2 / (
                            np.abs(np.log10(index[6] / index[7] + 1e-3)))
                    if OIII_5007_obs == 0 or SII_6725_obs == 0:
                        CHI_O3S2 = 0
                    elif index[6] == 0 or index[8] == 0:
                        CHI_O3S2 = tol_max
                    else:
                        CHI_O3S2 = (np.log10(index[6] / index[8]) - O3S2_obs) ** 2 / (
                            np.abs(np.log10(index[6] / index[8] + 1e-3)))

                    if ROIII_obs > 0:
                        CHI_OH = (CHI_ROIII ** 2 + CHI_NII ** 2 + CHI_OII ** 2 + CHI_OIII ** 2) ** 0.5
                    elif ROIII_obs == 0 and NII_6584_obs > 0:
                        CHI_OH = (CHI_NII ** 2 + CHI_O2O3 ** 2 + CHI_R23 ** 2 + CHI_O3N2 ** 2 + CHI_O3S2 ** 2) ** 0.5
                    elif ROIII_obs == 0 and NII_6584_obs == 0 and OIII_5007_obs > 0:
                        CHI_OH = (CHI_O2O3 ** 2 + CHI_R23 ** 2 + CHI_O3S2 ** 2) ** 0.5
                    elif ROIII_obs == 0 and OIII_5007_obs == 0:
                        CHI_OH = (CHI_O2Ne3 ** 2 + CHI_R2Ne3 ** 2) ** 0.5

                    if CHI_OH == 0:
                        OH_p = OH_p
                        logU_p = logU_p
                        den_OH = den_OH
                    else:
                        OH_p = index[0] / (CHI_OH) + OH_p
                        logU_p = index[2] / (CHI_OH) + logU_p
                        den_OH = 1 / (CHI_OH) + den_OH

                OH = OH_p / den_OH
                logU = logU_p / den_OH

            # Calculation of error of O/H and logU

            if R23_obs == -10 and NII_6584_obs == 0 and ROIII_obs == 0 and R2Ne3_obs == -10 and O3S2_obs == -10:
                eOH = 0
                elogU = 0
            else:
                CHI_ROIII = 0
                CHI_NII = 0
                CHI_OIII = 0
                CHI_OII = 0
                CHI_O2O3 = 0
                CHI_R23 = 0
                CHI_O2Ne3 = 0
                CHI_R2Ne3 = 0
                CHI_O3N2 = 0
                CHI_O3S2 = 0
                CHI_OH = 0
                for index in grid_c:
                    if ROIII_obs == 0:
                        CHI_ROIII = 0
                    elif index[5] == 0:
                        CHI_ROIII = tol_max
                    else:
                        CHI_ROIII = (index[6] / index[5] - ROIII_obs) ** 2 / (index[6] / index[5])
                    if OIII_5007_obs == 0:
                        CHI_OIII = 0
                    elif index[6] == 0:
                        CHI_OIII = tol_max
                    else:
                        CHI_OIII = (index[6] - OIII_5007_obs) ** 2 / index[6]
                    if OII_3727_obs == 0:
                        CHI_OII = 0
                    elif index[3] == 0:
                        CHI_OII = tol_max
                    else:
                        CHI_OII = (index[3] - OII_3727_obs) ** 2 / index[3]
                    if NII_6584_obs == 0:
                        CHI_NII = 0
                    elif index[7] == 0:
                        CHI_NII = tol_max
                    else:
                        CHI_NII = (index[7] - NII_6584_obs) ** 2 / index[7]
                    if OII_3727_obs == 0 or OIII_5007_obs == 0:
                        CHI_O2O3 = 0
                        CHI_R23 = 0
                    elif index[3] == 0 or index[6] == 0:
                        CHI_O2O3 = tol_max
                        CHI_R23 = tol_max
                    else:
                        CHI_O2O3 = (index[3] / index[6] - O2O3_obs) ** 2 / (index[3] / index[6])
                        CHI_R23 = (np.log10(index[3] + index[6]) - R23_obs) ** 2 / (
                            np.abs(np.log10(index[3] + index[6] + 1e-3)))
                    if OII_3727_obs == 0 or NeIII_3868_obs == 0:
                        CHI_O2Ne3 = 0
                        CHI_R2Ne3 = 0
                    elif index[3] == 0 or index[4] == 0:
                        CHI_O2Ne3 = tol_max
                        CHI_R2Ne3 = tol_max
                    else:
                        CHI_O2Ne3 = (index[3] / index[4] - O2Ne3_obs) ** 2 / (index[3] / index[4])
                        CHI_R2Ne3 = (np.log10(index[3] + index[4]) - R2Ne3_obs) ** 2 / (
                            np.abs(np.log10(index[3] + index[4] + 1e-3)))
                    if OIII_5007_obs == 0 or NII_6584_obs == 0:
                        CHI_O3N2 = 0
                    elif index[6] == 0 or index[7] == 0:
                        CHI_O3N2 = tol_max
                    else:
                        CHI_O3N2 = (np.log10(index[6] / index[7]) - O3N2_obs) ** 2 / (
                            np.abs(np.log10(index[6] / index[7] + 1e-3)))
                    if OIII_5007_obs == 0 or SII_6725_obs == 0:
                        CHI_O3S2 = 0
                    elif index[6] == 0 or index[8] == 0:
                        CHI_O3S2 = tol_max
                    else:
                        CHI_O3S2 = (np.log10(index[6] / index[8]) - O3S2_obs) ** 2 / (
                            np.abs(np.log10(index[6] / index[8] + 1e-3)))

                    if ROIII_obs > 0:
                        CHI_OH = (CHI_ROIII ** 2 + CHI_NII ** 2 + CHI_OII ** 2 + CHI_OIII ** 2) ** 0.5
                    elif ROIII_obs == 0 and NII_6584_obs > 0:
                        CHI_OH = (CHI_NII ** 2 + CHI_O2O3 ** 2 + CHI_R23 ** 2 + CHI_O3N2 ** 2 + CHI_O3S2 ** 2) ** 0.5
                    elif ROIII_obs == 0 and NII_6584_obs == 0 and OIII_5007_obs > 0:
                        CHI_OH = (CHI_O2O3 ** 2 + CHI_R23 ** 2 + CHI_O3S2 ** 2) ** 0.5
                    else:
                        CHI_OH = (CHI_O2Ne3 ** 2 + CHI_R2Ne3 ** 2) ** 0.5

                    if CHI_OH == 0:
                        OH_e = OH_e
                        logU_e = logU_e
                        den_OH_e = den_OH_e
                    else:
                        OH_e = (index[0] - OH) ** 2 / (CHI_OH) + OH_e
                        logU_e = (index[2] - logU) ** 2 / (CHI_OH) + logU_e
                        den_OH_e = 1 / (CHI_OH) + den_OH_e

                eOH = OH_e / den_OH_e
                elogU = logU_e / den_OH_e

            # Iterations for interpolated models
            if inter == 0 or OH == 0:
                OHf = OH
                logUf = logU
            elif inter == 1:
                igrid = interpolate_orig(grid_c, 2, logU - elogU - 0.25, logU + elogU + 0.25, 10)
                igrid = igrid[np.lexsort((igrid[:, 1], igrid[:, 2]))]
                igrid = interpolate_orig(igrid, 0, OH - eOH - 0.1, OH + eOH + 0.1, 10)
                igrid = igrid[np.lexsort((igrid[:, 0], igrid[:, 2]))]

                CHI_ROIII = 0
                CHI_OIII = 0
                CHI_OII = 0
                CHI_NII = 0
                CHI_O2O3 = 0
                CHI_R23 = 0
                CHI_O3N2 = 0
                CHI_O2Ne3 = 0
                CHI_R2Ne3 = 0
                CHI_O3S2 = 0
                CHI_OH = 0
                OH_p = 0
                logU_p = 0
                den_OH = 0

                for index in igrid:
                    if ROIII_obs == 0:
                        CHI_ROIII = 0
                    elif index[5] == 0:
                        CHI_ROIII = tol_max
                    else:
                        CHI_ROIII = (index[6] / index[5] - ROIII_obs) ** 2 / (index[6] / index[5])
                        if OIII_5007_obs == 0:
                            CHI_OIII = 0
                        elif index[6] == 0:
                            CHI_OIII = tol_max
                        else:
                            CHI_OIII = (index[6] - OIII_5007_obs) ** 2 / index[6]
                        if OII_3727_obs == 0:
                            CHI_OII = 0
                        elif index[3] == 0:
                            CHI_OII = tol_max
                        else:
                            CHI_OII = (index[3] - OII_3727_obs) ** 2 / index[3]
                    if NII_6584_obs == 0:
                        CHI_NII = 0
                    elif index[7] == 0:
                        CHI_NII = tol_max
                    else:
                        CHI_NII = (index[7] - NII_6584_obs) ** 2 / index[7]
                    if OII_3727_obs == 0 or OIII_5007_obs == 0:
                        CHI_O2O3 = 0
                        CHI_R23 = 0
                    elif index[3] == 0 or index[6] == 0:
                        CHI_O2O3 = tol_max
                        CHI_R23 = tol_max
                    else:
                        CHI_O2O3 = (index[3] / index[6] - O2O3_obs) ** 2 / (index[3] / index[6])
                        CHI_R23 = (np.log10(index[3] + index[6]) - R23_obs) ** 2 / (
                            np.abs(np.log10(index[3] + index[6] + 1e-3)))
                    if OII_3727_obs == 0 or NeIII_3868_obs == 0:
                        CHI_O2Ne3 = 0
                        CHI_R2Ne3 = 0
                    elif index[3] == 0 or index[4] == 0:
                        CHI_O2Ne3 = tol_max
                        CHI_R2Ne3 = tol_max
                    else:
                        CHI_O2Ne3 = (index[3] / index[4] - O2Ne3_obs) ** 2 / (index[3] / index[4])
                        CHI_R2Ne3 = (np.log10(index[3] + index[4]) - R2Ne3_obs) ** 2 / (
                            np.abs(np.log10(index[3] + index[4] + 1e-3)))
                    if OIII_5007_obs == 0 or NII_6584_obs == 0:
                        CHI_O3N2 = 0
                    elif index[6] == 0 or index[7] == 0:
                        CHI_O3N2 = tol_max
                    else:
                        CHI_O3N2 = (np.log10(index[6] / index[7]) - O3N2_obs) ** 2 / (
                            np.abs(np.log10(index[6] / index[7] + 1e-3)))
                    if OIII_5007_obs == 0 or SII_6725_obs == 0:
                        CHI_O3S2 = 0
                    elif index[6] == 0 or index[8] == 0:
                        CHI_O3S2 = tol_max
                    else:
                        CHI_O3S2 = (np.log10(index[6] / index[8]) - O3S2_obs) ** 2 / (
                            np.abs(np.log10(index[6] / index[8] + 1e-3)))

                    if ROIII_obs > 0:
                        CHI_OH = (CHI_ROIII ** 2 + CHI_NII ** 2 + CHI_OII ** 2 + CHI_OIII ** 2) ** 0.5
                    elif NII_6584_obs > 0 and OII_3727_obs > 0:
                        CHI_OH = (CHI_NII ** 2 + CHI_O2O3 ** 2 + CHI_R23 ** 2 + CHI_O2Ne3 ** 2 + CHI_R2Ne3 ** 2) ** 0.5
                    elif NII_6584_obs > 0 and OII_3727_obs == 0:
                        CHI_OH = (CHI_NII ** 2 + CHI_O3N2 ** 2 + CHI_O3S2 ** 2) ** 0.5
                    elif NII_6584_obs == 0:
                        CHI_OH = (CHI_O2O3 ** 2 + CHI_R23 ** 2 + CHI_O2Ne3 ** 2 + CHI_R2Ne3 ** 2 + CHI_O3S2 ** 2) ** 0.5

                    OH_p = index[0] / CHI_OH ** 2 + OH_p
                    logU_p = index[2] / CHI_OH ** 2 + logU_p
                    den_OH = 1 / CHI_OH ** 2 + den_OH

                if OH == 0:
                    OHf = OH
                    logUf = logU
                else:
                    OHf = OH_p / den_OH
                    logUf = logU_p / den_OH

            OH_mc.append(OHf)
            logU_mc.append(logUf)
            OHe_mc.append(eOH)
            logUe_mc.append(elogU)

        OHff = np.mean(OH_mc)
        eOHff = (np.std(OH_mc) ** 2 + np.mean(OHe_mc) ** 2) ** 0.5
        logUff = np.mean(logU_mc)
        elogUff = (np.std(logU_mc) ** 2 + np.std(logUe_mc) ** 2) ** 0.5

        OHffs.append(OHff)
        eOHffs.append(eOHff)
        NOffs.append(NOff)
        eNOffs.append(eNOff)
        logUffs.append(logUff)
        elogUffs.append(elogUff)

        if input0.size == 1 and tab == 0: continue
        print(round(100 * (count) / float(len(input1)), 1), '%', Names[tab], grid_type, '', round(OHff, 2),
              round(eOHff, 2),
              '', round(NOff, 2), round(eNOff, 2), '', round(logUff, 2), round(elogUff, 2))

    output['grid'] = grids
    output['OH'] = OHffs
    output['eOH'] = eOHffs
    output['NO'] = NOffs
    output['eNO'] = eNOffs
    output['logU'] = logUffs
    output['elogU'] = elogUffs

    output['grid'] = grids
    output['OH'] = OHffs
    output['eOH'] = eOHffs
    output['NO'] = NOffs
    output['eNO'] = eNOffs
    output['logU'] = logUffs
    output['elogU'] = elogUffs

    mc_dict = {'logU': np.array(logU_mc),
               'OH': np.array(OH_mc),
               'NO': np.array(NO_mc)}

    # if input0.size == 1:  output = np.delete(output, obj=1, axis=0)

    # lineas_header = [' HII-CHI-mistry v.5.1 output file', f'Input file: array',
    #                  'Iterations for MonteCarlo: ' + str(n),
    #                  'Used models: ' + sed_type, '',
    #                  'ID   O2Hb eO2Hb  Ne3Hb  eNeHb O3aHb  eO3aHb O3nHb  eO3nHb N2Hb   eN2Hb  S2Hb   eS2Hb  i O/H     eO/H  N/O    eN/O  logU   elogU']
    #
    # header = '\n'.join(lineas_header)
    #
    # np.savetxt(output_file, output, fmt=' '.join(['%s'] * 1 + ['%.3f'] * 12 + ['%i'] + ['%.2f'] * 6), header=header)
    # print('________________________________')
    # print('Results are stored in ' + str(output_file))

    return output