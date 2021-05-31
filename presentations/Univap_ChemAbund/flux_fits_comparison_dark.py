import numpy as np
import src.specsiser as sr
from src.specsiser.data_printing import PdfPrinter, latex_labels, numberStringFormat
from pathlib import Path
import pylatex

def percent_func(current_value, ref_value):
    return f' $({(1 - (current_value / ref_value)) * 100:.1f}\%)$'


conf_file_address = '/home/vital/PycharmProjects/vital_tests/astro/papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
resultsFolder = Path(obsData['file_information']['results_folder'])
tables_folder = Path('/home/vital/Dropbox/Astrophysics/Seminars/UniVapo 2021/')

ext = 'BR'
cycle = 'it3'
combined_line_dict = {'O2_3726A_m': 'O2_3726A-O2_3729A',
                      'O2_7319A_m': 'O2_7319A-O2_7330A',
                      'S2_6716A_m': 'S2_6716A-S2_6731A'}

# obj = 'gp121903'
for i, obj in enumerate(objList):

    if i < 3:

        objFolder = resultsFolder / f'{obj}'
        results_file = objFolder / f'{obj}_{ext}_measurements.txt'
        lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
        outputDb = objFolder / f'{obj}_{ext}_fitting_{cycle}.db'
        outputPiDb = objFolder / f'{obj}_{ext}_pI_fitting_{cycle}.db'
        outputDirectHIIchemDb = objFolder / f'{obj}_{ext}_Direct-Teff-logU_{cycle}.db'

        obj1_model = sr.SpectraSynthesizer()

        measurements_dict = sr.loadConfData(results_file, group_variables=False)
        objLinesDF = sr.import_emission_line_data(lineLog_file)
        pI_results = sr.load_MC_fitting(outputPiDb)
        fit_results = sr.load_MC_fitting(outputDb)
        bayes_plus_PI_results = sr.load_MC_fitting(outputDirectHIIchemDb)

        # Recover Direct method Bayesian fitting
        sr_lineLabels = fit_results['Input_data']['lineLabels_list']
        inFlux = fit_results['Input_data']['inputFlux_array']
        inErr = fit_results['Input_data']['inputErr_array']
        flux_matrix = fit_results['trace']['calcFluxes_Op']
        mean_line_values = flux_matrix.mean(axis=0)
        std_line_values = flux_matrix.std(axis=0)

        # Recover Direct method + HIIchem fitting
        srPlus_lineLabels = bayes_plus_PI_results['Input_data']['lineLabels_list']
        inFluxPlus = bayes_plus_PI_results['Input_data']['inputFlux_array']
        inErrPlus = bayes_plus_PI_results['Input_data']['inputErr_array']
        flux_matrix_plus = bayes_plus_PI_results['trace']['calcFluxes_Op']
        mean_plus_line_values = flux_matrix_plus.mean(axis=0)
        std_plus_line_values = flux_matrix_plus.std(axis=0)

        # Recover Photo-ionization Bayesian fitting
        pI_lineLabels = pI_results['Input_data']['lineLabels_list']
        inFlux_PI = pI_results['Input_data']['inputFlux_array']
        inErr_PI = pI_results['Input_data']['inputErr_array']
        # flux_matrix_PI = pI_results['trace']['calcFluxes_Op']
        # mean_linePi_values = flux_matrix_PI.mean(axis=0)
        # std_linePi_values = flux_matrix_PI.std(axis=0)
        flux_tensor_dict = {}
        trace = pI_results['trace']
        for i_line, line in enumerate(pI_lineLabels):
            tensor_ref = f'{pI_lineLabels[i_line]}_Op'
            if tensor_ref in trace.varnames:
                if line == 'S2_6716A':
                    flux_tensor_dict['S2_6716A_m'] = trace[tensor_ref]
                else:
                    flux_tensor_dict[line] = trace[tensor_ref]

        # Add missing rows
        S2_6716A_m_flux = objLinesDF.loc['S2_6716A', 'obsFlux'] + objLinesDF.loc['S2_6731A', 'obsFlux']
        S2_6716A_m_err = np.sqrt(objLinesDF.loc['S2_6716A', 'obsFluxErr']**2 + objLinesDF.loc['S2_6731A', 'obsFluxErr']**2)
        S2_6716A_m_lambda = (objLinesDF.loc['S2_6716A', 'obsFlux'] + objLinesDF.loc['S2_6731A', 'obsFlux']) / 2
        objLinesDF.loc['S2_6716A_m', 'obsFlux'] = S2_6716A_m_flux
        objLinesDF.loc['S2_6716A_m', 'obsFluxErr'] = S2_6716A_m_err
        objLinesDF.loc['S2_6716A_m', 'f_lambda'] = S2_6716A_m_lambda

        tableLines = sr_lineLabels.copy()
        i_6725A = sr_lineLabels.index('S2_6731A') + 1
        tableLines.insert(i_6725A, 'S2_6716A_m')

        ion_array, wave_array, latexLabel_array = sr.label_decomposition(tableLines, combined_dict=combined_line_dict)

        table_headers = ['Line ID',
                         'Observed flux',
                         r'\makecell{Bayesian: \\ Direct method}',
                         r'\makecell{\textsc{HII-CHI-mistry: } \\ Spherical}',
                         r'\makecell{\textsc{HII-CHI-mistry: } \\ Plane-Parallel}',
                         r'\makecell{Bayesian \textsc{HII-CHI-mistry: } \\ Plane-Parallel}',
                         r'\makecell{Bayesian \\ Direct Method + \textsc{HII-CHI-mistry} \\ Plane-Parallel}'
                         ]

        pdf = PdfPrinter()
        pdf.create_pdfDoc(pdf_type='table')
        pdf.pdfDoc.append(pylatex.NoEscape('\definecolor{background}{rgb}{0.169, 0.169, 0.169}'))
        pdf.pdfDoc.append(pylatex.NoEscape('\definecolor{foreground}{rgb}{0.702, 0.780, 0.847}'))
        pdf.pdf_insert_table(table_headers, addfinalLine=True, color_background='background', color_font='foreground')

        for i, lineLabel in enumerate(tableLines):

            row_data = ['-'] * len(table_headers)
            row_data[0] = latexLabel_array[i]

            obsFlux, obsErr = objLinesDF.loc[lineLabel, 'obsFlux'], objLinesDF.loc[lineLabel, 'obsFluxErr']
            row_data[1] = r'${:.3f}\pm{:.3f}$'.format(obsFlux, obsErr)

            # Bayesian direct method fitting
            if lineLabel in sr_lineLabels:
                i_line = sr_lineLabels.index(lineLabel)
                # row_data[1] =  r'${:.3f}\pm{:.3f}$'.format(inFlux[i_line], inErr[i_line])
                row_data[2] = r'${:.3f}\pm{:.3f}$'.format(mean_line_values[i_line], std_line_values[i_line])
                row_data[2] += percent_func(obsFlux, mean_line_values[i_line])

            # Spherical fitting HII-CHI-mistry
            cHbeta_label = obsData[obj]['cHbeta_label']
            cHbeta = np.array(measurements_dict[f'Extinction_{cycle}'][cHbeta_label], dtype=float)
            flambda = objLinesDF.loc[lineLabel, 'f_lambda']
            print(lineLabel, obsFlux, obsFlux * np.power(10, cHbeta[0] * flambda), obsErr * np.power(10, cHbeta[0] * flambda))

            section = f'HII_Tef_fit_C17_bb_Teff_30-90_sph.dat_{cycle}_sed'
            if lineLabel in measurements_dict[section]:
                lineInt, lineIntErr = np.array(measurements_dict[section][lineLabel], dtype=float)
                lineFlux, lineFluxErr = lineInt / np.power(10, cHbeta[0] * flambda), lineIntErr / np.power(10, cHbeta[0] * flambda)
                row_data[3] = r'${:.3f}\pm{:.3f}$'.format(lineFlux, lineFluxErr)
                row_data[3] += percent_func(obsFlux, lineFlux)

            section = f'HII_Tef_fit_C17_bb_Teff_30-90_pp.dat_{cycle}_sed'
            if lineLabel in measurements_dict[section]:
                lineInt, lineIntErr = np.array(measurements_dict[section][lineLabel], dtype=float)
                lineFlux, lineFluxErr = lineInt / np.power(10, cHbeta[0] * flambda), lineIntErr / np.power(10, cHbeta[0] * flambda)
                row_data[4] = r'${:.3f}\pm{:.3f}$'.format(lineFlux, lineFluxErr)
                row_data[4] += percent_func(obsFlux, lineFlux)

            if lineLabel in flux_tensor_dict:
                lineFlux, lineFluxErr = np.mean(flux_tensor_dict[lineLabel]), np.std(flux_tensor_dict[lineLabel])
                row_data[5] = r'${:.3f}\pm{:.3f}$'.format(lineFlux, lineFluxErr)
                row_data[5] += percent_func(obsFlux, lineFlux)

            if lineLabel in srPlus_lineLabels:
                i_line = sr_lineLabels.index(lineLabel)
                row_data[6] = r'${:.3f}\pm{:.3f}$'.format(mean_plus_line_values[i_line], std_plus_line_values[i_line])
                row_data[6] += percent_func(obsFlux, mean_plus_line_values[i_line])

            pdf.addTableRow(row_data, last_row=False, color_background='background', color_font='foreground')

        pdf.table.add_hline()

        table_file = tables_folder/f'{obj}_lineFlux_fit_comparison'
        pdf.generate_pdf(table_file, clean_tex=False)


