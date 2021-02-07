import numpy as np
import pyneb as pn
import pandas as pd
import src.specsiser as sr
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt, rcParams
from src.specsiser.data_printing import label_decomposition, PdfPrinter
from lmfit.models import LinearModel
from lmfit import fit_report
from uncertainties import umath, unumpy, ufloat


FLUX_TEX_TABLE_HEADERS = [r'$Transition$', '$f(\lambda)$', '$F(\lambda)$', '$I(\lambda)$']
FLUX_TXT_TABLE_HEADERS = [r'Transition', 'f(lambda)', 'F(lambda)', 'F(lambda)_error', 'I(lambda)', 'I(lambda)_error']


def red_corr_HalphaHbeta_ratio(lines_df, default_cHbeta):

    # Normalizing flux
    if 'H1_6563A' in lines_df.index:
        flux_Halpha = lines_df.loc['H1_6563A', 'gauss_flux']
        flux_Hbeta = lines_df.loc['H1_4861A', 'intg_flux']
        halpha_norm = flux_Halpha / flux_Hbeta

        rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
        rc.setCorr(obs_over_theo=halpha_norm / 2.86, wave1=6563., wave2=4861.)
        cHbeta = rc.cHbeta
    else:
        flux_Hbeta = lines_df.loc['H1_4861A', 'intg_flux']
        rc = pn.RedCorr(R_V=3.4, law='G03 LMC', cHbeta=default_cHbeta)
        cHbeta = float(rc.cHbeta)

    return cHbeta, rc


def compute_cHbeta(input_labels, line_df, reddening_curve, R_v, temp=10000.0, den=100.0, ref_wave='H1_4861A'):

    assert ref_wave in line_df.index, f'- ERROR: {ref_wave} not found in input lines log dataframe for c(Hbeta) calculation'

    # Label the lines which are found in the lines log
    idcs_lines = line_df.index.isin(input_labels)
    line_labels = line_df.loc[idcs_lines].index.values
    ion_ref, waves_ref, latexLabels_ref = label_decomposition(ref_wave, scalar_output=True)
    ion_array, waves_array, latexLabels_array = label_decomposition(line_labels)

    # Observed ratios
    Href_flux, Href_err = line_df.loc[ref_wave, 'intg_flux'], line_df.loc[ref_wave, 'intg_err']
    obsFlux, obsErr = np.empty(line_labels.size), np.empty(line_labels.size)
    slice_df = line_df.loc[idcs_lines]
    idcs_intg = slice_df.blended == 'None'
    idcs_gauss = ~idcs_intg
    obsFlux[idcs_intg], obsErr[idcs_intg] = slice_df.loc[idcs_intg, 'intg_flux'], slice_df.loc[idcs_intg, 'intg_err']
    obsFlux[idcs_gauss], obsErr[idcs_gauss] = slice_df.loc[idcs_gauss, 'gauss_flux'], slice_df.loc[idcs_gauss, 'gauss_err']
    obsRatio_uarray = unumpy.uarray(obsFlux, obsErr)/ufloat(Href_flux, Href_err)

    # Theoretical ratios
    H1 = pn.RecAtom('H', 1)
    refEmis = H1.getEmissivity(tem=temp, den=den, wave=waves_ref)
    emisIterable = (H1.getEmissivity(tem=temp, den=den, wave=wave) for wave in waves_array)
    linesEmis = np.fromiter(emisIterable, float)
    theoRatios = linesEmis/refEmis

    # Reddening law
    rc = pn.RedCorr(R_V=R_v, law=reddening_curve)
    X_Hbeta, Xx_ref, Xx = rc.X(4861.0), rc.X(waves_ref), rc.X(waves_array)
    f_lines = Xx/Xx_ref - 1
    f_Hbeta = X_Hbeta/Xx_ref - 1

    # cHbeta linear fit values
    x_values = f_lines - f_Hbeta
    y_values = np.log10(theoRatios) - unumpy.log10(obsRatio_uarray)

    # Perform fit
    lineModel = LinearModel()
    y_nom, y_std = unumpy.nominal_values(y_values), unumpy.std_devs(y_values)
    pars = lineModel.make_params(intercept=y_nom.min(), slope=0)
    output = lineModel.fit(y_nom, pars, x=x_values, weights=1/np.sqrt(y_std))
    cHbeta, cHbeta_err = output.params['slope'].value, output.params['slope'].stderr
    intercept, intercept_err = output.params['intercept'].value, output.params['intercept'].stderr

    output_dict = dict(cHbeta=cHbeta, cHbeta_err=cHbeta_err,
                       intercept=intercept, intercept_err=intercept_err,
                       obsRecomb=unumpy.nominal_values(obsRatio_uarray), obsRecombErr=unumpy.std_devs(obsRatio_uarray),
                       y=y_nom, y_err=y_std,
                       x=x_values)

    return output_dict


def compute_arms_flambda(line_DF, red_law, R_V, ref_line='H1_4861A'):

    required_columns = {'wavelength'}
    assert required_columns.issubset(line_DF.columns), f'-- ERROR: one or severals columns from {required_columns} missing in lines log'
    assert ref_line in line_DF.index, f'- ERROR: {ref_line} not found in line log dataframe for f_lambda calculation'

    wave_ref = line_DF.loc[ref_line].wavelength
    wavelengths = line_DF.wavelength

    # Reddening law
    rc = pn.RedCorr(R_V=R_V, law=red_law)
    X_ref, Xx = rc.X(wave_ref), rc.X(wavelengths)
    f_ref, f_x = X_ref/X_ref - 1, Xx/X_ref - 1

    return f_ref, f_x


def compute_spectrum_flambda(wavelength, red_law, R_V, ref_line='H1_4861A'):

    # TODO combine previous two functions to work wiht reddening corrections for both continuum and lines
    if ref_line == 'H1_4861A':
        wave_ref = 4861.0
    elif ref_line == 'H1_6563A':
        wave_ref = 6563.0

    # Reddening law
    rc = pn.RedCorr(R_V=R_V, law=red_law)
    X_ref, Xx = rc.X(wave_ref), rc.X(wavelength)
    f_ref, f_x = X_ref/X_ref - 1, Xx/X_ref - 1

    return f_x, f_ref


def double_arm_redCorr(wave_spec, flux_spec, wave_boundary, red_law, red_R_V, cHbeta):

    rc = pn.RedCorr(R_V=red_R_V, law=red_law, cHbeta=cHbeta[0])

    # Mix correction
    idcs_waveBlue = wave_spec < wave_boundary
    idcs_waveRed = wave_spec > wave_boundary
    f_spec_blue, f_Hbeta = compute_spectrum_flambda(wave_spec[idcs_waveBlue], rc.law, rc.R_V, ref_line='H1_4861A')
    f_spec_red, f_Halpha = compute_spectrum_flambda(wave_spec[idcs_waveRed],  rc.law, rc.R_V, ref_line='H1_6563A')

    corr_array = np.zeros(wave_spec.size)
    corr_array[idcs_waveBlue] = np.power(10, cHbeta[0] * f_spec_blue) * np.power(10, cHbeta[0])
    corr_array[idcs_waveRed] = np.power(10, cHbeta[0] * f_spec_red) * np.power(10, 0.4 * rc.E_BV * rc.X(6563))

    int_array = flux_spec * corr_array

    return int_array, corr_array


def normalize_flux(line_DF, norm_line, scale_factor=1, flux_mode='auto'):

    required_columns = {'blended', 'intg_flux', 'gauss_flux', 'intg_err', 'gauss_err', 'latexLabel'}
    assert required_columns.issubset(line_DF.columns), f'-- ERROR: one or severals columns from {required_columns} missing in lines log'
    assert norm_line in line_DF.index, f'\n-- ERROR: normalizing line {norm_line} was not found in lines log'

    # Normalization wave
    flux_label = 'intg_' if line_DF.loc[norm_line, 'blended'] == 'None' else 'gauss_'
    norm_flux, norm_err = line_DF.loc[norm_line, [f'{flux_label}flux', f'{flux_label}err']]

    # Assign wavelengths # TODO add strong method to distinguish merged lines whose desired flux is the integrated
    idcs_gauss = (line_DF.blended != 'None') & (~line_DF.latexLabel.str.contains('+', regex=False))
    idcs_ingt = ~idcs_gauss
    obsFlux, obsErr = np.zeros(len(line_DF.index)), np.zeros(len(line_DF.index))
    obsFlux[idcs_gauss], obsFlux[idcs_ingt] = line_DF.loc[idcs_gauss, 'gauss_flux'], line_DF.loc[idcs_ingt, 'intg_flux']
    obsErr[idcs_gauss], obsErr[idcs_ingt] = line_DF.loc[idcs_gauss, 'gauss_err'], line_DF.loc[idcs_ingt, 'intg_err']

    # Scale the fluxes
    obsFlux, obsErr = obsFlux * scale_factor, obsErr * scale_factor

    # Normalize fluxes
    obsFlux_uarray = unumpy.uarray(obsFlux, obsErr) / ufloat(norm_flux, norm_err)

    return unumpy.nominal_values(obsFlux_uarray), unumpy.std_devs(obsFlux_uarray)


def deredd_fluxes(obs_flux, obs_err, cHbeta_nom, cHbeta_err, lines_flambda):

    # Generate uncertainty variables to propagate the error
    cHbeta = ufloat(cHbeta_nom, cHbeta_err),
    obsFlux_uarray = unumpy.uarray(obs_flux, obs_err)

    # Compute line intensity
    obsInt_uarray = obsFlux_uarray * unumpy.pow(10, cHbeta * lines_flambda)
    obsInt, obsIntErr = unumpy.nominal_values(obsInt_uarray), unumpy.std_devs(obsInt_uarray)

    return obsInt, obsIntErr


def table_fluxes(lineLabels, f_lambda, flux, flux_err, inten, inten_err, cHbeta, cHbeta_err, ref_label, ref_flux, ref_err,
                 output_address, scaleTable=1000):

    # TODO this could be included in sr.print
    txt_address = f'{output_address}.txt'
    pdf_address = f'{output_address}'

    # Measure line fluxes
    pdf = PdfPrinter()
    pdf.create_pdfDoc(pdf_address, pdf_type='table')
    pdf.pdf_insert_table(FLUX_TEX_TABLE_HEADERS)

    # Dataframe as container as a txt file
    tableDF = pd.DataFrame(columns=FLUX_TXT_TABLE_HEADERS[1:])

    # Get the reference
    ion_array, wavelength_array, latexLabel_array = label_decomposition(lineLabels)
    ion_ref, wavelength_ref, latexLabel_ref = label_decomposition(ref_label, scalar_output=True)

    for i, linelabel in enumerate(lineLabels):

        flambda_entry = f'{f_lambda[i]:0.2f}'
        flux_entry = r'${:0.2f}\,\pm\,{:0.2f}$'.format(flux[i] * scaleTable, flux_err[i] * scaleTable)
        intensity_entry = r'${:0.2f}\,\pm\,{:0.2f}$'.format(inten[i] * scaleTable, inten_err[i] * scaleTable)

        # Add row of data
        pdf_row_i = [latexLabel_array[i], flambda_entry, flux_entry, intensity_entry]
        txt_row_i = [linelabel, f_lambda[i],flux[i], flux_err[i], inten[i], inten_err[i]]

        lastRow_check = True if linelabel == lineLabels[-1] else False
        pdf.addTableRow(pdf_row_i, last_row=lastRow_check)
        tableDF.loc[linelabel] = txt_row_i[1:]

    # Data last rows
    row_Hbetaflux = [r'$H\beta$ $(erg\,cm^{-2} s^{-1} \AA^{-1})$',
                     '',
                     r'${:0.3e}\,\pm\,{:0.3e}$'.format(ref_flux, ref_err),
                     '']

    row_cHbeta = [r'$c(H\beta)$',
                  '',
                  r'${:0.2f}\,\pm\,{:0.2f}$'.format(cHbeta, cHbeta_err),
                  '']

    pdf.addTableRow(row_Hbetaflux, last_row=False)
    pdf.addTableRow(row_cHbeta, last_row=False)
    tableDF.loc[row_Hbetaflux[0]] = row_Hbetaflux[1:] + [''] * (tableDF.columns.size - 3)
    tableDF.loc[row_cHbeta[0]] = row_cHbeta[1:] + [''] * (tableDF.columns.size - 3)

    # Format last rows
    pdf.table.add_hline()
    pdf.table.add_hline()

    # Save the pdf table
    try:
        pdf.generate_pdf(clean_tex=True)
    except:
        print('-- PDF compilation failure')

    # Save the txt table
    with open(txt_address, 'wb') as output_file:
        string_DF = tableDF.to_string()
        string_DF = string_DF.replace('$', '')
        output_file.write(string_DF.encode('UTF-8'))

    return


def exitinction_corr_plot(objName, corr_dict_list, ext_file_list, plot_save_file=None):

    ext_color_dict = dict(BR='tab:purple', B='tab:blue', SDSS='black')

    STANDARD_PLOT = {'figure.figsize': (14, 7),
                     'axes.titlesize': 14,
                     'axes.labelsize': 18,
                     'legend.fontsize': 12,
                     'xtick.labelsize': 12,
                     'ytick.labelsize': 12}
    rcParams.update(STANDARD_PLOT)

    fig, ax = plt.subplots(figsize=(8, 4))

    ax.update({'xlabel': r'$f_{\lambda} - f_{H\beta}$',
               'ylabel': r'$ \left(\frac{I_{\lambda}}{I_{\H\beta}}\right)_{Theo} - \left(\frac{F_{\lambda}}{F_{\H\beta}}\right)_{Obs}$',
               'title': f' {objName} logaritmic extinction calculation'})

    for idx_file, corr_lines_dict in enumerate(corr_dict_list):

        ext_file = ext_file_list[idx_file]
        corr_dict, blue_rc = corr_lines_dict['four_lines'], corr_lines_dict['three_lines']

        ax.errorbar(corr_dict['x'], corr_dict['y'], yerr=corr_dict['y_err'], color=ext_color_dict[ext_file], fmt='o')
        all_ylineFit = corr_dict['cHbeta'] * corr_dict['x'] + corr_dict['intercept']
        label = r'$c(H\beta)$ = ${:.2f}\pm{:.2f}$ ' \
                r'($H\alpha$, $H\beta$, $H\gamma$, $H\delta$)'.format(corr_dict['cHbeta'], corr_dict['cHbeta_err'])
        ax.plot(corr_dict['x'], all_ylineFit, color=ext_color_dict[ext_file], label=label, linestyle='--')


        blue_ylineFit = blue_rc['cHbeta'] * blue_rc['x'] + blue_rc['intercept']
        label = r'$c(H\beta)$ = ${:.2f}\pm{:.2f}$ ($H\beta$, $H\gamma$, $H\delta$)'.format(blue_rc['cHbeta'],
                                                                                           blue_rc['cHbeta_err'])
        ax.plot(blue_rc['x'], blue_ylineFit, color=ext_color_dict[ext_file], label=label, linestyle=':')

    ax.legend()
    plt.tight_layout()
    if plot_save_file is not None:
        plt.savefig(plot_save_file, dpi=200, bbox_inches='tight')
    else:
        plt.show()

    return