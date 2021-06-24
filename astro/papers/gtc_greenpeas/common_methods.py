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


def compute_cHbeta(lines_dict, line_df, reddening_curve, R_v, temp=10000.0, den=100.0, ref_wavesList=['H1_4861A'],
                   compMode='auto', plot_address=None):

    ref_wavesList = np.array(ref_wavesList, ndmin=1)

    H1 = pn.RecAtom('H', 1)

    output_dict = {}

    for i, item in enumerate(lines_dict.items()):

        cHbeta_label, input_labels = item
        ref_wave = ref_wavesList[i] if ref_wavesList.size > 1 else ref_wavesList[0]
        assert ref_wave in line_df.index, f'- ERROR: {ref_wave} not found in input lines log dataframe for c(Hbeta) calculation'

        # Label the lines which are found in the lines log
        idcs_lines = line_df.index.isin(input_labels)
        line_labels = line_df.loc[idcs_lines].index.values
        ion_ref, waves_ref, latexLabels_ref = label_decomposition(ref_wave, scalar_output=True)
        ion_array, waves_array, latexLabels_array = label_decomposition(line_labels)

        # Observed ratios
        if compMode == 'auto':
            Href_flux, Href_err = line_df.loc[ref_wave, 'intg_flux'], line_df.loc[ref_wave, 'intg_err']
            obsFlux, obsErr = np.empty(line_labels.size), np.empty(line_labels.size)
            slice_df = line_df.loc[idcs_lines]
            idcs_intg = slice_df.blended == 'None'
            # obsFlux[idcs_intg], obsErr[idcs_intg] = slice_df.loc[idcs_intg, ['intg_flux', 'intg_err']].values
            # obsFlux[~idcs_intg], obsErr[~idcs_intg] = slice_df.loc[~idcs_intg, ['gauss_flux', 'gauss_err']]
            obsFlux[idcs_intg] = slice_df.loc[idcs_intg, 'intg_flux'].values
            obsErr[idcs_intg] = slice_df.loc[idcs_intg, 'intg_err'].values
            obsFlux[~idcs_intg] = slice_df.loc[~idcs_intg, 'gauss_flux'].values
            obsErr[~idcs_intg] = slice_df.loc[~idcs_intg, 'gauss_err'].values
            obsRatio_uarray = unumpy.uarray(obsFlux, obsErr) / ufloat(Href_flux, Href_err)

        elif compMode == 'gauss':
            Href_flux, Href_err = line_df.loc[ref_wave, 'gauss_flux'], line_df.loc[ref_wave, 'gauss_err']
            obsFlux, obsErr = line_df.loc[idcs_lines, 'gauss_flux'], line_df.loc[idcs_lines, 'gauss_err']
            obsRatio_uarray = unumpy.uarray(obsFlux, obsErr) / ufloat(Href_flux, Href_err)

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

        # cHbeta linear fit values
        x_values = f_lines - f_ref
        y_values = np.log10(theoRatios) - unumpy.log10(obsRatio_uarray)

        # Perform fit
        lineModel = LinearModel()
        y_nom, y_std = unumpy.nominal_values(y_values), unumpy.std_devs(y_values)
        pars = lineModel.make_params(intercept=y_nom.min(), slope=0)
        output = lineModel.fit(y_nom, pars, x=x_values, weights=1 / np.sqrt(y_std))
        cHbeta, cHbeta_err = output.params['slope'].value, output.params['slope'].stderr
        intercept, intercept_err = output.params['intercept'].value, output.params['intercept'].stderr

        # Store the results
        output_dict[cHbeta_label] = dict(cHbeta=cHbeta,
                                         cHbeta_err=cHbeta_err,
                                         intercept=intercept,
                                         intercept_err=intercept_err,
                                         obsRecomb=unumpy.nominal_values(obsRatio_uarray),
                                         obsRecombErr=unumpy.std_devs(obsRatio_uarray),
                                         y=y_nom, y_err=y_std,
                                         x=x_values)

    if plot_address is not None:

        STANDARD_PLOT = {'figure.figsize': (14, 7), 'axes.titlesize': 12, 'axes.labelsize': 14,
                         'legend.fontsize': 10, 'xtick.labelsize': 10, 'ytick.labelsize': 10}
        rcParams.update(STANDARD_PLOT)

        axes_dict = {'xlabel': r'$f_{\lambda} - f_{H\beta}$',
                     'ylabel': r'$ \left(\frac{I_{\lambda}}{I_{\H\beta}}\right)_{Theo} - \left(\frac{F_{\lambda}}{F_{\H\beta}}\right)_{Obs}$',
                     'title': f'Logaritmic extinction coefficient calculation'}

        fig, ax = plt.subplots(figsize=(8, 4))
        fig.subplots_adjust(bottom=-0.7)

        for i, item in enumerate(lines_dict.items()):

            cHbeta_label, input_labels = item
            ext_fit_results = output_dict[cHbeta_label]

            ion_list, waves_list, latexLabels_list = label_decomposition(input_labels, scalar_output=True)

            data_points = ax.errorbar(ext_fit_results['x'], ext_fit_results['y'], yerr=ext_fit_results['y_err'],
                        fmt='o', color='black')

            label = r'{}: $c(H\beta)$ = ${:.2f}\pm{:.2f}$'.format(cHbeta_label, ext_fit_results['cHbeta'], ext_fit_results['cHbeta_err'])
            cHbeta_trndline = ext_fit_results['cHbeta'] * ext_fit_results['x'] + ext_fit_results['intercept']
            ax.plot(ext_fit_results['x'], cHbeta_trndline, linestyle='--', label=label)

        ax.update(axes_dict)
        # ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.15))

        # Shrink current axis's height by 10% on the bottom
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.7])

        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.25),
                  fancybox=True, shadow=True)

        plt.tight_layout()
        plt.savefig(plot_address, dpi=200, bbox_inches='tight')
        # plt.show()

    return output_dict


def compute_arms_flambda(line_DF, red_law, R_V, ref_line='H1_4861A'):
    required_columns = {'wavelength'}
    assert required_columns.issubset(
        line_DF.columns), f'-- ERROR: one or severals columns from {required_columns} missing in lines log'
    assert ref_line in line_DF.index, f'- ERROR: {ref_line} not found in line log dataframe for f_lambda calculation'

    wave_ref = line_DF.loc[ref_line].wavelength
    wavelengths = line_DF.wavelength

    # Reddening law
    rc = pn.RedCorr(R_V=R_V, law=red_law)
    X_ref, Xx = rc.X(wave_ref), rc.X(wavelengths)
    f_ref, f_x = X_ref / X_ref - 1, Xx / X_ref - 1

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
    f_ref, f_x = X_ref / X_ref - 1, Xx / X_ref - 1

    return f_x, f_ref


def double_arm_redCorr(wave_spec, flux_spec, wave_boundary, red_law, red_R_V, cHbeta):
    rc = pn.RedCorr(R_V=red_R_V, law=red_law, cHbeta=cHbeta[0])

    # Mix correction
    idcs_waveBlue = wave_spec < wave_boundary
    idcs_waveRed = wave_spec > wave_boundary
    f_spec_blue, f_Hbeta = compute_spectrum_flambda(wave_spec[idcs_waveBlue], rc.law, rc.R_V, ref_line='H1_4861A')
    f_spec_red, f_Halpha = compute_spectrum_flambda(wave_spec[idcs_waveRed], rc.law, rc.R_V, ref_line='H1_6563A')

    corr_array = np.zeros(wave_spec.size)
    corr_array[idcs_waveBlue] = np.power(10, cHbeta[0] * f_spec_blue) * np.power(10, cHbeta[0])
    corr_array[idcs_waveRed] = np.power(10, cHbeta[0] * f_spec_red) * np.power(10, 0.4 * rc.E_BV * rc.X(6563))

    int_array = flux_spec * corr_array

    return int_array, corr_array


def normalize_flux(line_DF, norm_line, scale_factor=1, flux_mode='auto'):
    required_columns = {'blended_label', 'intg_flux', 'gauss_flux', 'intg_err', 'gauss_err', 'latexLabel'}
    assert required_columns.issubset(
        line_DF.columns), f'-- ERROR: one or severals columns from {required_columns} missing in lines log'
    assert norm_line in line_DF.index, f'\n-- ERROR: normalizing line {norm_line} was not found in lines log'

    # Normalization wave
    flux_label = 'intg_' if line_DF.loc[norm_line, 'blended_label'] == 'None' else 'gauss_'
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


def table_fluxes(lineLabels, f_lambda, flux, flux_err, inten, inten_err, cHbeta, cHbeta_err, ref_label, ref_flux,
                 ref_err,
                 output_address, scaleTable=1000):
    # TODO this could be included in sr.print
    txt_address = f'{output_address}.txt'
    pdf_address = f'{output_address}'

    # Measure line fluxes
    pdf = PdfPrinter()
    pdf.create_pdfDoc(pdf_type='table')
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
        txt_row_i = [linelabel, f_lambda[i], flux[i], flux_err[i], inten[i], inten_err[i]]

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
        pdf.generate_pdf(pdf_address, clean_tex=True)
    except:
        print('-- PDF compilation failure')

    # Save the txt table
    with open(txt_address, 'wb') as output_file:
        string_DF = tableDF.to_string()
        string_DF = string_DF.replace('$', '')
        output_file.write(string_DF.encode('UTF-8'))

    return


def exitinction_corr_plot(objName, corr_dict_list, ext_file_list, plot_save_file=None):
    STANDARD_PLOT = {'figure.figsize': (14, 7), 'axes.titlesize': 14, 'axes.labelsize': 18,
                     'legend.fontsize': 12, 'xtick.labelsize': 12, 'ytick.labelsize': 12}
    rcParams.update(STANDARD_PLOT)

    axes_dict = {'xlabel': r'$f_{\lambda} - f_{H\beta}$',
                 'ylabel': r'$ \left(\frac{I_{\lambda}}{I_{\H\beta}}\right)_{Theo} - \left(\frac{F_{\lambda}}{F_{\H\beta}}\right)_{Obs}$',
                 'title': f' {objName} logaritmic extinction calculation'}

    ext_color_dict = dict(BR='tab:purple', B='tab:blue', SDSS='black')

    fig, ax = plt.subplots(figsize=(8, 4))

    for iCor, corr_lines_dict in enumerate(corr_dict_list):
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

    ax.update(axes_dict)
    ax.legend()
    plt.tight_layout()
    if plot_save_file is not None:
        plt.savefig(plot_save_file, dpi=200, bbox_inches='tight')
    else:
        plt.show()

    return


def exitinction_corr_plot_backUp(objName, corr_dict_list, ext_file_list, plot_save_file=None):
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


def interpolate(grid, z, zmin, zmax, n):

   ncol = 10
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
         elif grid[row,x] < grid[row+1,x] or grid[row,y] < grid[row+1,y] :
            vec.append(grid[row,col])
            no_inter = no_inter + 1
         else:
            inter = inter + 1
            for index in range(0,n):
               i = grid[row,col]+(index)*(grid[row+1,col]-grid[row,col])/n
               vec.append(i)
   out = np.transpose(np.reshape(vec,(-1,n*inter+no_inter)))
   return out


def epm_fitting(input00, output_file, n, sed, geo, inter, grid_file=None, HCm_folder=''):

    input0 = np.genfromtxt(input00, dtype=None, names=True)

    if input0.size == 1:
        input1 = np.stack((input0, input0))
    else:
        input1 = input0

    if geo == 1 and sed == 1:
        bin = 99
        fileAddress = f'{HCm_folder}C17_WMb_Teff_30-60_pp.dat'
        grid = np.loadtxt(fileAddress)
        if inter == 0:
            sed_type = 'WM-Basic stellar atmosphere. Plane-parallel geometry. Not interpolated'
            print('Teff and U calculation using WM-Basic models with plane-paralell geometry and non-interpolation')
        elif inter == 1:
            sed_type = 'WM-Basic stellar atmosphere. Plane-parallel geometry. Interpolated'
            print('Teff and U calculation using WM-Basic models with plane-paralell geometry and interpolation')

    elif geo == 2 and sed == 1:
        bin = 99
        fileAddress = f'{HCm_folder}C17_WMb_Teff_30-60_sph.dat'
        grid = np.loadtxt(fileAddress)
        if inter == 0:
            sed_type = 'WM-Basic stellar atmosphereSpherical geometry. Not interpolated'
            print('Teff and U calculation using WM-Basic models with spherical geometry and non-interpolation')
        elif inter == 1:
            sed_type = 'WM-Basic stellar atmosphere. Spherical geometry. Interpolated'
            print('Teff and U calculation using WM-Basic models with spherical geometry and interpolation')

    elif geo == 1 and sed == 2:
        bin = 132
        fileAddress = f'{HCm_folder}C17_bb_Teff_30-90_pp.dat'
        grid = np.loadtxt(fileAddress)
        if inter == 0:
            sed_type = 'Black body. Plane-parallel geometry. Not interpolated'
            print('Teff and U calculation using black body models with plane-parallel geometry and non-interpolation')
        elif inter == 1:
            sed_type = 'Black body. Plane-parallel geometry. Interpolated'
            print('Teff and U calculation using black body models with plane-parallel geometry and interpolation')

    elif geo == 2 and sed == 2:
        bin = 132
        fileAddress = f'{HCm_folder}C17_bb_Teff_30-90_sph.dat'
        grid = np.loadtxt(fileAddress)
        if inter == 0:
            sed_type = 'Black body. spherical geometry. Not interpolated'
            print('Teff and U calculation using black body models with spherical geometry and non-interpolation')
        elif inter == 1:
            sed_type = 'Black body. spherical geometry. Interpolated'
            print('Teff and U calculation using black body models with spehrical geometry and interpolation')

    elif geo == 1 and sed == 3:
        sed_type = 'BPASS cluster atmospheres, Mup = 300. Plane-parallel geometry'
        bin = 63
        fileAddress = f'{HCm_folder}C17_bpass21_imf135_300_pp_esc.dat'
        grid = np.loadtxt(fileAddress)
        print('f_abs and U calculation using BPASS models with plane-parallel geometry')

    elif geo == 2 and sed == 3:
        bin = 77
        fileAddress = f'{HCm_folder}C17_bpass21_imf135_300_sph_esc.dat'
        grid = np.loadtxt(fileAddress)
        if inter == 0:
            sed_type = 'BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. Spherical geometry. Not interpolated'
            print('F_abs and U calculation usingBPASS models with spherical geometry and non-interpolation')
        elif inter == 1:
            sed_type = 'BPASS cluster atmospheres, Mup = 300, x = 1.35, age = 4 Myr. Spherical geometry. Interpolated'
            print('F_abs and U calculation usingBPASS models with spherical geometry and interpolation')

    if grid_file is not None:
        grid = np.loadtxt(grid_file)

    print(f'Using grids: {grid_file}')

    OHffs = []
    eOHffs = []
    Teffs = []
    eTeffs = []
    logUffs = []
    elogUffs = []

    Label_ID = False
    Label_OH = False
    Label_eOH = False
    Label_OII = False
    Label_eOII = False
    Label_OIII_4959 = False
    Label_eOIII_4959 = False
    Label_OIII_5007 = False
    Label_eOIII_5007 = False
    Label_SII = False
    Label_eSII = False
    Label_SII_6716 = False
    Label_eSII_6716 = False
    Label_SII_6731 = False
    Label_eSII_6731 = False
    Label_SIII_9069 = False
    Label_eSIII_9069 = False
    Label_SIII_9532 = False
    Label_eSIII_9532 = False
    Label_HeI_4471 = False
    Label_eHeI_4471 = False
    Label_HeI_5876 = False
    Label_eHeI_5876 = False
    Label_HeII_4686 = False
    Label_eHeII_4686 = False

    for col in range(0, len(input1.dtype.names), 1):
        if input1.dtype.names[col] == 'ID':
            Label_ID = True
        if input1.dtype.names[col] == '12logOH':
            Label_OH = True
        if input1.dtype.names[col] == 'e12logOH':
            Label_eOH = True
        if input1.dtype.names[col] == 'OII_3727':
            Label_OII = True
        if input1.dtype.names[col] == 'eOII_3727':
            Label_eOII = True
        if input1.dtype.names[col] == 'OIII_4959':
            Label_OIII_4959 = True
        if input1.dtype.names[col] == 'eOIII_4959':
            Label_eOIII_4959 = True
        if input1.dtype.names[col] == 'OIII_5007':
            Label_OIII_5007 = True
        if input1.dtype.names[col] == 'eOIII_5007':
            Label_eOIII_5007 = True
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
        if input1.dtype.names[col] == 'SIII_9069':
            Label_SIII_9069 = True
        if input1.dtype.names[col] == 'eSIII_9069':
            Label_eSIII_9069 = True
        if input1.dtype.names[col] == 'SIII_9532':
            Label_SIII_9532 = True
        if input1.dtype.names[col] == 'eSIII_9532':
            Label_eSIII_9532 = True
        if input1.dtype.names[col] == 'HeI_4471':
            Label_HeI_4471 = True
        if input1.dtype.names[col] == 'eHeI_4471':
            Label_eHeI_4471 = True
        if input1.dtype.names[col] == 'HeI_5876':
            Label_HeI_5876 = True
        if input1.dtype.names[col] == 'eHeI_5876':
            Label_eHeI_5876 = True
        if input1.dtype.names[col] == 'HeII_4686':
            Label_HeII_4686 = True
        if input1.dtype.names[col] == 'eHeII_4686':
            Label_eHeII_4686 = True

    if Label_ID == False:
        Names = np.arange(1, input1.size + 1, 1)
    else:
        Names = input1['ID']
    if Label_OH == False:
        logOH = np.zeros(input1.size)
    else:
        logOH = input1['12logOH']
    if Label_eOH == False:
        elogOH = np.zeros(input1.size)
    else:
        elogOH = input1['e12logOH']
    if Label_OII == False:
        OII_3727 = np.zeros(input1.size)
    else:
        OII_3727 = input1['OII_3727']
    if Label_eOII == False:
        eOII_3727 = np.zeros(input1.size)
    else:
        eOII_3727 = input1['eOII_3727']
    if Label_OIII_4959 == False and Label_OIII_5007 == False:
        OIII_5007 = np.zeros(input1.size)
    elif Label_OIII_4959 == False and Label_OIII_5007 == True:
        OIII_5007 = input1['OIII_5007']
    elif Label_OIII_4959 == True and Label_OIII_5007 == False:
        OIII_5007 = 4 * input1['OIII_4959']
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
    if Label_SIII_9069 == False and Label_SIII_9532 == False:
        SIII_9069 = np.zeros(input1.size)
    elif Label_SIII_9069 == False and Label_SIII_9532 == True:
        SIII_9069 = input1['SIII_9069'] / 2.44
    elif Label_SIII_9069 == True and Label_SIII_9532 == False:
        SIII_9069 = input1['SIII_9069']
    else:
        SIII_9069 = (input1['SIII_9069'] + input1['SIII_9532']) / 3.44
    if Label_eSIII_9069 == False and Label_eSIII_9532 == False:
        eSIII_9069 = np.zeros(input1.size)
    elif Label_eSIII_9069 == False and Label_eSIII_9532 == True:
        eSIII_9069 = input1['eSIII_9069'] / 2.44
    elif Label_eSIII_9069 == True and Label_eSIII_9532 == False:
        eSIII_9069 = input1['eSIII_9069']
    else:
        eSIII_9069 = (input1['eSIII_9069'] + input1['eSIII_9532']) / 3.44
    if Label_HeI_4471 == False:
        HeI_4471 = np.zeros(input1.size)
    else:
        HeI_4471 = input1['HeI_4471']
    if Label_eHeI_4471 == False:
        eHeI_4471 = np.zeros(input1.size)
    else:
        eHeI_4471 = input1['eHeI_4471']
    if Label_HeI_5876 == False:
        HeI_5876 = np.zeros(input1.size)
    else:
        HeI_5876 = input1['HeI_5876']
    if Label_eHeI_5876 == False:
        eHeI_5876 = np.zeros(input1.size)
    else:
        eHeI_5876 = input1['eHeI_5876']
    if Label_HeII_4686 == False:
        HeII_4686 = np.zeros(input1.size)
    else:
        HeII_4686 = input1['HeII_4686']
    if Label_eHeII_4686 == False:
        eHeII_4686 = np.zeros(input1.size)
    else:
        eHeII_4686 = input1['eHeII_4686']

    output = np.zeros(input1.size,
                      dtype=[('ID', 'U12'), ('OII_3727', float), ('eOII_3727', float), ('OIII_5007', float),
                             ('eOIII_5007', float), ('HeI_4471', float), ('eHeI_4471', float), ('HeI_5876', float),
                             ('eHeI_5876', float), ('HeII_4686', float), ('eHeII_4686', float), ('SII_6725', float),
                             ('eSII_6725', float), ('SIII_9069', float), ('eSIII_9069', float), ('OH', float),
                             ('eOH', float), ('Teff', float), ('eTeff', float), ('logU', float), ('elogU', float)])

    output['ID'] = Names
    output['OII_3727'] = OII_3727
    output['eOII_3727'] = eOII_3727
    output['OIII_5007'] = OIII_5007
    output['eOIII_5007'] = eOIII_5007
    output['HeI_4471'] = HeI_4471
    output['eHeI_4471'] = eHeI_4471
    output['HeI_5876'] = HeI_5876
    output['eHeI_5876'] = eHeI_5876
    output['HeII_4686'] = HeII_4686
    output['eHeII_4686'] = eHeII_4686
    output['SII_6725'] = SII_6725
    output['eSII_6725'] = eSII_6725
    output['SIII_9069'] = SIII_9069
    output['eSIII_9069'] = eSIII_9069

    print('Reading grids ....')
    print('')
    print('')
    print('---------------------------------------')
    if sed < 3:
        print('(%)   ID.   12+lyog(O/H)  T_eff(K)    log(U)')
    else:
        print('(%)   ID     12+log(O/H)  f_abs       log(U)')

    print('---------------------------------------')

    # Beginning of loop of calculation

    count = 0
    for tab in range(0, len(input1), 1):

        count = count + 1

        OH_mc = []
        Teff_mc = []
        logU_mc = []
        eOH_mc = []
        eTeff_mc = []
        elogU_mc = []

        for monte in range(0, n, 1):

            OH_p = 0
            logU_p = 0
            Teff_p = 0
            den_OH = 0
            den_Teff = 0
            OH_e = 0
            Teff_e = 0
            logU_e = 0
            den_OH_e = 0
            den_Teff_e = 0
            tol_max = 1e2

            OII_3727_obs = 0
            if OII_3727[tab] > 0:
                while OII_3727_obs <= 0:
                    OII_3727_obs = np.random.normal(OII_3727[tab], eOII_3727[tab] + 1e-3)
            OIII_5007_obs = 0
            if OIII_5007[tab] > 0:
                while OIII_5007_obs <= 0:
                    OIII_5007_obs = np.random.normal(OIII_5007[tab], eOIII_5007[tab] + 1e-3)
            SII_6725_obs = 0
            if SII_6725[tab] > 0:
                while SII_6725_obs <= 0:
                    SII_6725_obs = np.random.normal(SII_6725[tab], eSII_6725[tab] + 1e-3)
            SIII_9069_obs = 0
            if SIII_9069[tab] > 0:
                while SIII_9069_obs <= 0:
                    SIII_9069_obs = np.random.normal(SIII_9069[tab], eSIII_9069[tab] + 1e-3)
            HeI_4471_obs = 0
            if HeI_4471[tab] > 0:
                while HeI_4471_obs <= 0:
                    HeI_4471_obs = np.random.normal(HeI_4471[tab], eHeI_4471[tab] + 1e-3)
            HeI_5876_obs = 0
            if HeI_5876[tab] > 0:
                while HeI_5876_obs <= 0:
                    HeI_5876_obs = np.random.normal(HeI_5876[tab], eHeI_5876[tab] + 1e-3)
            HeII_4686_obs = 0
            if HeII_4686[tab] > 0:
                while HeII_4686_obs <= 0:
                    HeII_4686_obs = np.random.normal(HeII_4686[tab], eHeII_4686[tab] + 1e-3)
            if OII_3727_obs == 0 or OIII_5007_obs == 0:
                O2O3_obs = -10
                R23_obs = -10
            else:
                O2O3_obs = np.log10(OII_3727_obs / OIII_5007_obs)
                R23_obs = np.log10(OII_3727_obs + OIII_5007_obs)

            if SII_6725_obs == 0 or SIII_9069_obs == 0:
                S2S3_obs = -10
                S23_obs = -10
            else:
                S2S3_obs = np.log10(SII_6725_obs / SIII_9069_obs)
                S23_obs = (SII_6725_obs + SIII_9069_obs)
            if SII_6725_obs == 0 or OIII_5007_obs == 0:
                S2O3_obs = -10
            else:
                S2O3_obs = np.log10(SII_6725_obs / OIII_5007_obs)
            if HeI_4471_obs == 0 or HeII_4686_obs == 0:
                He12a_obs = -10
            else:
                He12a_obs = np.log10(HeI_4471_obs / HeII_4686_obs)
            if HeI_5876_obs == 0 or HeII_4686_obs == 0:
                He12b_obs = -10
            else:
                He12b_obs = np.log10(HeI_5876_obs / HeII_4686_obs)

            # Interpolation of grid at specific O/H

            if logOH[tab] > 0:
                OH = np.random.normal(logOH[tab], elogOH[tab] + 1e-3)
                OH_mc.append(OH)

                grid_T0 = []
                if OH <= 7.1:
                    OH = 7.1
                    i0 = 0
                    i1 = bin
                elif OH >= 7.1 and OH < 7.4:
                    i0 = 0
                    i1 = bin
                elif OH >= 7.4 and OH < 7.7:
                    i0 = bin
                    i1 = 2 * bin
                elif OH >= 7.7 and OH < 8.0:
                    i0 = 2 * bin
                    i1 = 3 * bin
                elif OH >= 8.0 and OH < 8.3:
                    i0 = 3 * bin
                    i1 = 4 * bin
                elif OH >= 8.3 and OH < 8.6:
                    i0 = 4 * bin
                    i1 = 5 * bin
                elif OH >= 8.6 and OH < 8.9:
                    i0 = 5 * bin
                    i1 = 6 * bin
                elif OH >= 8.9:
                    OH = 8.9
                    i0 = 5 * bin
                    i1 = 6 * bin

                for x in range(0, bin):
                    for y in range(0, 10):
                        grid_T0.append(
                            grid[i0 + x, y] * np.abs(0.3 - OH + grid[i0, 0]) / 0.3 + grid[i1 + x, y] * np.abs(
                                0.3 - grid[i1, 0] + OH) / 0.3)

                #         grid_T0.append(grid[i0+x,y]*np.abs(0.3-grid[i0,0]+OH)/0.3 + grid[i1+x,y]*np.abs(0.3-grid[i1,0]+OH)/0.3)

                grid_T = np.reshape(grid_T0, (bin, 10))

            else:
                OH = 0
                OH_mc.append(OH)
                grid_T = grid

            np.savetxt('int_models.dat', grid_T, fmt='%.2f')

            # Calculation of T and log U

            if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10:
                Teff = 0
                logU = 0
            else:
                CHI_O2O3 = 0
                CHI_S2S3 = 0
                CHI_S23 = 0
                CHI_S2O3 = 0
                CHI_He12a = 0
                CHI_He12b = 0
                CHI_He12 = 0

                for index in grid_T:
                    if index[9] == 0 and HeII_4686_obs > 0: continue
                    if S2S3_obs == -10:
                        CHI_S2S3 = 0
                        CHI_S23 = 0
                    elif index[5] == 0 or index[6] == 0:
                        CHI_S2S3 = tol_max
                        CHI_S23 = tol_max
                    else:
                        CHI_S2S3 = (np.log10(index[5] / index[6]) - S2S3_obs) ** 2 / S2S3_obs
                        CHI_S23 = (index[5] + index[6] - S23_obs) ** 2 / S23_obs
                    if S2O3_obs == -10:
                        CHI_S2O3 = 0
                    elif index[5] == 0 or index[4] == 0:
                        CHI_S2O3 = tol_max
                    else:
                        CHI_S2O3 = (np.log10(index[5] / index[4]) - S2O3_obs) ** 2 / S2O3_obs
                    if O2O3_obs == -10:
                        CHI_O2O3 = 0
                    elif index[3] == 0 or index[4] == 0:
                        CHI_O2O3 = tol_max
                    else:
                        CHI_O2O3 = (np.log10(index[3] / index[4]) - O2O3_obs) ** 2 / O2O3_obs
                    if He12a_obs == -10:
                        CHI_He12a = 0
                    elif index[7] == 0 or index[9] == 0:
                        CHI_He12a = tol_max
                    else:
                        CHI_He12a = (np.log10(index[7] / index[9]) - He12a_obs) ** 2 / He12a_obs
                    if He12b_obs == -10:
                        CHI_He12b = 0
                    elif index[8] == 0 or index[9] == 0:
                        CHI_He12b = tol_max
                    else:
                        CHI_He12b = (np.log10(index[8] / index[9]) - He12b_obs) ** 2 / He12b_obs
                    if CHI_He12a == 0 and CHI_He12b == 0:
                        CHI_HE12 = 0
                    elif CHI_He12a == 0:
                        CHI_He12 = CHI_He12b
                    elif CHI_He12b == 0:
                        CHI_He12 = CHI_He12a
                    else:
                        CHI_He12 = (CHI_He12a + CHI_He12b) / 2

                    if OII_3727_obs == 0:
                        CHI_Teff = (CHI_S2S3 ** 2 + CHI_He12 ** 2 + CHI_S2O3 ** 2) ** 0.5
                    elif SIII_9069_obs == 0 and HeII_4686_obs == 0:
                        CHI_Teff = (CHI_O2O3 ** 2 + CHI_S2O3 ** 2) ** 0.5
                    elif SIII_9069_obs == 0 and HeII_4686_obs > 0:
                        CHI_Teff = (CHI_He12 ** 2 + CHI_O2O3 ** 2) ** 0.5 # Estuve aqui
                    else:
                        CHI_Teff = (CHI_S2S3 ** 2 + CHI_O2O3 ** 2 + CHI_He12 ** 2) ** 0.5

                    Teff_p = index[1] * (1 / CHI_Teff) ** 2 + Teff_p
                    logU_p = index[2] * (1 / CHI_Teff) ** 2 + logU_p
                    den_Teff = (1 / CHI_Teff) ** 2 + den_Teff
                Teff = Teff_p / den_Teff
                logU = logU_p / den_Teff

            # Calculation of T and log U errors

            if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10:
                eTeff = 0
                elogU = 0
            else:
                CHI_O2O3 = 0
                CHI_S2S3 = 0
                CHI_S23 = 0
                CHI_S2O3 = 0
                CHI_He12a = 0
                CHI_He12b = 0
                CHI_He12 = 0

                for index in grid_T:
                    if index[9] == 0 and HeII_4686_obs > 0: continue
                    if S2S3_obs == -10:
                        CHI_S2S3 = 0
                        CHI_S23 = 0
                    elif index[5] == 0 or index[6] == 0:
                        CHI_S2S3 = tol_max
                        CHI_S23 = tol_max
                    else:
                        CHI_S2S3 = (np.log10(index[5] / index[6]) - S2S3_obs) ** 2 / S2S3_obs
                        CHI_S23 = (index[5] + index[6] - S23_obs) ** 2 / S23_obs
                    if S2O3_obs == -10:
                        CHI_S2O3 = 0
                    elif index[5] == 0 or index[4] == 0:
                        CHI_S2O3 = tol_max
                    else:
                        CHI_S2O3 = (np.log10(index[5] / index[4]) - S2O3_obs) ** 2 / S2O3_obs
                    if O2O3_obs == -10:
                        CHI_O2O3 = 0
                    elif index[3] == 0 or index[4] == 0:
                        CHI_O2O3 = tol_max
                    else:
                        CHI_O2O3 = (np.log10(index[3] / index[4]) - O2O3_obs) ** 2 / O2O3_obs
                    if He12a_obs == -10:
                        CHI_He12a = 0
                    elif index[7] == 0 or index[9] == 0:
                        CHI_He12a = tol_max
                    else:
                        CHI_He12a = (np.log10(index[7] / index[9]) - He12a_obs) ** 2 / He12a_obs
                    if He12b_obs == -10:
                        CHI_He12b = 0
                    elif index[8] == 0 or index[9] == 0:
                        CHI_He12b = tol_max
                    else:
                        CHI_He12b = (np.log10(index[8] / index[9]) - He12b_obs) ** 2 / He12b_obs
                    if CHI_He12a == 0 and CHI_He12b == 0:
                        CHI_HE12 = 0
                    elif CHI_He12a == 0:
                        CHI_He12 = CHI_He12b
                    elif CHI_He12b == 0:
                        CHI_He12 = CHI_He12a
                    else:
                        CHI_He12 = (CHI_He12a + CHI_He12b) / 2

                    if OII_3727_obs == 0:
                        CHI_Teff = (CHI_S2S3 ** 2 + CHI_He12 ** 2 + CHI_S2O3 ** 2) ** 0.5
                    elif SIII_9069_obs == 0 and HeII_4686_obs == 0:
                        CHI_Teff = (CHI_O2O3 ** 2 + CHI_S2O3 ** 2) ** 0.5
                    elif SIII_9069_obs == 0 and HeII_4686_obs > 0:
                        CHI_Teff = (CHI_He12 ** 2 + CHI_O2O3 ** 2) ** 0.5
                    else:
                        CHI_Teff = (CHI_S2S3 ** 2 + CHI_O2O3 ** 2 + CHI_He12 ** 2) ** 0.5

                    if sed < 3:
                        Teff_e = np.abs(index[1] - Teff) * (1 / CHI_Teff) ** 2 + Teff_e
                    else:
                        Teff_e = np.abs(np.log10(index[1] + 1e-5) - np.log10(Teff)) * (1 / CHI_Teff) ** 2 + Teff_e
                    logU_e = np.abs(index[2] - logU) * (1 / CHI_Teff) ** 2 + logU_e
                    den_Teff_e = 1 * (1 / CHI_Teff) ** 2 + den_Teff_e

                eTeff = Teff_e / den_Teff_e
                if sed == 3:
                    eTeff = Teff * np.log10(eTeff)
                elogU = logU_e / den_Teff_e

                # Iterations for the interpolation mode

                if inter == 0:
                    Teff = Teff
                    logU = logU
                elif inter == 1:
                    igrid = grid_T[np.lexsort((grid_T[:, 0], grid_T[:, 2]))]
                    igrid = interpolate(igrid, 1, Teff - eTeff - 1e4, Teff + eTeff + 1e4, 10)
                    igrid = igrid[np.lexsort((igrid[:, 0], igrid[:, 1]))]
                    igrid = interpolate(igrid, 2, logU - elogU - 0.25, logU + elogU + 0.25, 10)

                    #            np.savetxt('int_models.dat',igrid,fmt='%.2f')

                    if S2S3_obs == -10 and O2O3_obs == -10 and He12a_obs == -10 and He12b_obs == -10:
                        Teff = 0
                        logU = 0
                    else:
                        CHI_O2O3 = 0
                        CHI_S2S3 = 0
                        CHI_S23 = 0
                        CHI_S2O3 = 0
                        CHI_He12a = 0
                        CHI_He12b = 0
                        CHI_He12 = 0

                        for index in igrid:
                            if index[9] == 0 and HeII_4686_obs > 0: continue
                            if S2S3_obs == -10:
                                CHI_S2S3 = 0
                                CHI_S23 = 0
                            elif index[5] == 0 or index[6] == 0:
                                CHI_S2S3 = tol_max
                                CHI_S23 = tol_max
                            else:
                                CHI_S2S3 = (np.log10(index[5] / index[6]) - S2S3_obs) ** 2 / S2S3_obs
                                CHI_S23 = (index[5] + index[6] - S23_obs) ** 2 / S23_obs
                            if S2O3_obs == -10:
                                CHI_S2O3 = 0
                            elif index[5] == 0 or index[4] == 0:
                                CHI_S2O3 = tol_max
                            else:
                                CHI_S2O3 = (np.log10(index[5] / index[4]) - S2O3_obs) ** 2 / S2O3_obs
                            if O2O3_obs == -10:
                                CHI_O2O3 = 0
                            elif index[3] == 0 or index[4] == 0:
                                CHI_O2O3 = tol_max
                            else:
                                CHI_O2O3 = (np.log10(index[3] / index[4]) - O2O3_obs) ** 2 / O2O3_obs
                            if He12a_obs == -10:
                                CHI_He12a = 0
                            elif index[7] == 0 or index[9] == 0:
                                CHI_He12a = tol_max
                            else:
                                CHI_He12a = (np.log10(index[7] / index[9]) - He12a_obs) ** 2 / He12a_obs
                            if He12b_obs == -10:
                                CHI_He12b = 0
                            elif index[8] == 0 or index[9] == 0:
                                CHI_He12b = tol_max
                            else:
                                CHI_He12b = (np.log10(index[8] / index[9]) - He12b_obs) ** 2 / He12b_obs
                            if CHI_He12a == 0 and CHI_He12b == 0:
                                CHI_HE12 = 0
                            elif CHI_He12a == 0:
                                CHI_He12 = CHI_He12b
                            elif CHI_He12b == 0:
                                CHI_He12 = CHI_He12a
                            else:
                                CHI_He12 = (CHI_He12a + CHI_He12b) / 2

                            if OII_3727_obs == 0:
                                CHI_Teff = (CHI_S2S3 ** 2 + CHI_He12 ** 2 + CHI_S2O3 ** 2) ** 0.5
                            elif SIII_9069_obs == 0 and HeII_4686_obs == 0:
                                CHI_Teff = (CHI_O2O3 ** 2 + CHI_S2O3 ** 2) ** 0.5
                            elif SIII_9069_obs == 0 and HeII_4686_obs > 0:
                                CHI_Teff = (CHI_He12 ** 2 + CHI_O2O3 ** 2) ** 0.5
                            else:
                                CHI_Teff = (CHI_S2S3 ** 2 + CHI_O2O3 ** 2 + CHI_He12 ** 2) ** 0.5

                            if sed < 3:
                                Teff_p = index[1] * (1 / CHI_Teff) + Teff_p
                            else:
                                Teff_p = index[1] * (1 / CHI_Teff) + Teff_p
                            logU_p = index[2] * (1 / CHI_Teff) + logU_p
                            den_Teff = (1 / CHI_Teff) + den_Teff

                        Teff = Teff_p / den_Teff
                        if sed == 3:
                            Teff = Teff
                        logU = logU_p / den_Teff

            Teff_mc.append(Teff)
            logU_mc.append(logU)
            eTeff_mc.append(eTeff)
            elogU_mc.append(elogU)

        if logOH[tab] > 0:
            OHf = logOH[tab]
            eOHf = elogOH[tab]
        else:
            OHf = 0
            eOHf = 0
        Tefff = np.mean(Teff_mc)
        eTefff = (np.std(Teff_mc) ** 2 + np.mean(eTeff_mc) ** 2) ** 0.5
        logUf = np.mean(logU_mc)
        elogUf = (np.std(logU_mc) ** 2 + np.mean(elogU_mc) ** 2) ** 0.5

        OHffs.append(OHf)
        eOHffs.append(eOHf)
        Teffs.append(Tefff)
        eTeffs.append(eTefff)
        logUffs.append(logUf)
        elogUffs.append(elogUf)

        if input0.size == 1 and tab == 0: continue

        if sed == 3:
            print(round(100 * (count) / float(len(input1)), 1), '%', '', Names[tab], round(OHf, 2), round(eOHf, 2),
                  round(Tefff, 2), round(eTefff, 2), round(logUf, 2), round(elogUf, 2))
        else:
            print(round(100 * (count) / float(len(input1)), 1), '%', '', Names[tab], round(OHf, 2), round(eOHf, 2),
                  100 * int(Tefff / 100), 100 * int(eTefff / 100), round(logUf, 2), round(elogUf, 2))

    output['OH'] = OHffs
    output['eOH'] = eOHffs
    output['Teff'] = Teffs
    output['eTeff'] = eTeffs
    output['logU'] = logUffs
    output['elogU'] = elogUffs

    if input0.size == 1:  output = np.delete(output, obj=1, axis=0)

    if sed == 3:
        lineas_header = [' HII-CHI-mistry v.5.0 output file', ' Input file:' + input00,
                         'Iterations for MonteCarlo: ' + str(n), 'Used models: ' + sed_type, '',
                         'ID. O2Hb eO2Hb O3Hb  eO3Hb  eO3Hb S2Hb  eS2Hb S3Hb  eS3Hb 4471Hb e4471Hb 5678Hb e5648Hb He2Hb eHe2Hb O/H   eO/H  fabs  efabs logU elogU']
    else:
        lineas_header = [' HII-CHI-mistry-Teff v.5.0 output file', ' Input file:' + input00,
                         'Iterations for MonteCarlo: ' + str(n), 'Used models: ' + sed_type, '',
                         'ID.   O2Hb eO2Hb O3Hb  eO3Hb 4471Hb e4471Hb 5876Hb e5876Hb He2Hb eHe2Hb  S2Hb  eS2Hb S3Hb  eS3Hb  O/H   eO/H  Teff  eTeff logU elogU']

    header = '\n'.join(lineas_header)

    if sed == 3:
        np.savetxt(output_file, output,
                   fmt=' '.join(['%s'] * 1 + ['%.3f'] * 14 + ['%.2f'] * 2 + ['%.2f'] * 2 + ['%.2f'] * 2), header=header)
    else:
        np.savetxt(output_file, output,
                   fmt=' '.join(['%s'] * 1 + ['%.3f'] * 14 + ['%.2f'] * 2 + ['%.0f'] * 2 + ['%.2f'] * 2), header=header)
    print('________________________________')
    print('Results are stored in ' + input00 + '_hcm-teff-output.dat')

    return


def reading_epm_grids():

    grid_file = 'D:/Dropbox/Astrophysics/Tools/HCm-Teff_v5.01/C17_bb_Teff_30-90_pp.dat'

    lineConversionDict = dict(O2_3726A_m='OII_3727',
                           O3_5007A='OIII_5007',
                           He1_4471A='HeI_4471',
                           He1_5876A='HeI_5876',
                           He2_4686A='HeII_4686',
                           S2_6716A='SII_6717,31',
                           S3_9069A='SIII_9069')

    # Load the data and get axes range
    grid_array = np.loadtxt(grid_file)

    grid_axes = dict(OH=np.unique(grid_array[:, 0]),
                     Teff=np.unique(grid_array[:, 1]),
                     logU=np.unique(grid_array[:, 2]))

    # Sort the array according to 'logU', 'Teff', 'OH'
    idcs_sorted_grid = np.lexsort((grid_array[:, 1], grid_array[:, 2], grid_array[:, 0]))
    sorted_grid = grid_array[idcs_sorted_grid]

    # Loop throught the emission line and abundances and restore the grid
    grid_dict = {}
    for i, item in enumerate(lineConversionDict.items()):
        lineLabel, epmLabel = item

        grid_dict[lineLabel] = np.zeros((grid_axes['logU'].size,
                                         grid_axes['Teff'].size,
                                         grid_axes['OH'].size))

        for j, abund in enumerate(grid_axes['OH']):

            idcsSubGrid = sorted_grid[:, 0] == abund
            lineGrid = sorted_grid[idcsSubGrid, i + 3]
            lineMatrix = lineGrid.reshape((grid_axes['logU'].size, grid_axes['Teff'].size))
            grid_dict[lineLabel][:, :, j] = lineMatrix[:, :]

    return grid_dict, grid_axes


def check_previous_measurements(objName, parameter_list, measurements_dict, cfg_dict):

    output_dict = {}

    # Chemical abundances
    for param in parameter_list:
        if param in measurements_dict['it3_ionic_Abundances']:
            output_dict[param] = measurements_dict['it3_ionic_Abundances'][param]

    if 'O2_3726A_m' in measurements_dict['it3_ionic_Abundances']:
        output_dict['O2'] = measurements_dict['it3_ionic_Abundances']['O2_3726A_m']

    conversion_dict = {'He1': 'He1r', 'He2': 'He2r'}
    for fit_label, measure_label in conversion_dict.items():
        if measure_label in measurements_dict['it3_ionic_Abundances']:
                output_dict[fit_label] = measurements_dict['it3_ionic_Abundances'][measure_label]

    # Extinction
    cHbeta_key = cfg_dict[objName]['cHbeta_label']
    output_dict['cHbeta'] = measurements_dict['Extinction_it3'][cHbeta_key]

    # Electron parameters
    conversion_dict = {'n_e': 'ne', 'T_low': 'Te_low', 'T_high': 'Te_high'}
    for fit_label, measure_label in conversion_dict.items():
        if measure_label in measurements_dict['it3_electron_parameters']:
            output_dict[fit_label] = measurements_dict['it3_electron_parameters'][measure_label]

    return output_dict