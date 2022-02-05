import numpy as np
import pandas as pd
import configparser
import pyneb as pn

from pathlib import Path
from astropy.io import fits
from mpdaf.obj import Cube
from astropy.table import Table
from lime.tools import label_decomposition
from lmfit.models import LinearModel
from lime.io import format_for_table
from collections import Sequence

# State target lines and parameters
target_lines = ['H1_4861A', 'H1_4861A_w1', 'H1_6563A',  'H1_6563A_w1', 'H1_6563A_w2', 'H1_8750A', 'H1_8863A', 'H1_9015A', 'H1_9229A',
                'O3_4959A', 'O3_5007A', 'O3_5007A_w1',
                'S2_6716A', 'S2_6731A',
                'S3_6312A', 'S3_9069A',
                'He1_5876A', 'He1_6678A']

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
                'sigma_vel_err': target_lines}

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

grid_HII_CHI_mistry_conversion = {'logOH': '12+log(O/H)',
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


def voxel_security_check(linesDF):

    check = False

    if 'H1_4861A' in linesDF.index:
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