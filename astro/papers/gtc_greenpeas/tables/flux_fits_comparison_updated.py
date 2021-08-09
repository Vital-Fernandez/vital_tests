import numpy as np
import src.specsiser as sr
from src.specsiser.data_printing import PdfPrinter, latex_labels, numberStringFormat
from pathlib import Path
import pandas as pd

def percent_func(current_value, ref_value):
    return f' $({(1 - (current_value / ref_value)) * 100:.1f}\%)$'

def percent_func_txt(current_value, ref_value):
    return f'{(1 - (current_value / ref_value)) * 100:.1f}'

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
resultsFolder = Path(obsData['file_information']['results_folder'])
tables_folder = Path(obsData['file_information']['tables_folder'])

ext = 'BR'
cycle = 'it3'
combined_line_dict = {'O2_3726A_m': 'O2_3726A-O2_3729A',
                      'O2_7319A_m': 'O2_7319A-O2_7330A',
                      'S2_6716A_m': 'S2_6716A-S2_6731A'}

for i, obj in enumerate(objList):

    if i < 3:

        objFolder = resultsFolder / f'{obj}'
        results_file = objFolder / f'{obj}_{ext}_measurements.txt'
        lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
        outputDb = objFolder / f'{obj}_{ext}_fitting_{cycle}.db'
        outputTxt = objFolder / f'{obj}_{ext}_fitting_{cycle}.txt'
        outputPiDb = objFolder / f'{obj}_{ext}_pI_fitting_{cycle}.db'
        outputPiTxt = objFolder / f'{obj}_{ext}_pI_fitting_{cycle}.txt'
        outputDirectHIIchemDb = objFolder / f'{obj}_{ext}_Direct-Teff-logU_{cycle}.db'
        outputDirectHIIchemTxt = objFolder / f'{obj}_{ext}_Direct-Teff-logU_{cycle}.txt'

        obj1_model = sr.SpectraSynthesizer()

        measurements_dict = sr.loadConfData(results_file, group_variables=False)
        objLinesDF = sr.import_emission_line_data(lineLog_file)

        fit_pickle = sr.load_MC_fitting(outputDb)
        pI_pickle = sr.load_MC_fitting(outputPiDb)
        bayes_plus_PI_pickle = sr.load_MC_fitting(outputDirectHIIchemDb)

        fit_inputs_dict = sr.loadConfData(outputTxt)
        pI_inputs_dict = sr.loadConfData(outputPiTxt)
        bayes_inputs_dict = sr.loadConfData(outputDirectHIIchemTxt)

        # Recover Direct method Bayesian fitting
        sr_lineLabels = fit_inputs_dict['Input_data']['lineLabels_list']
        inFlux = fit_inputs_dict['Input_data']['inputFlux_array']
        inErr = fit_inputs_dict['Input_data']['inputErr_array']
        flux_matrix = fit_pickle['trace']['calcFluxes_Op']
        mean_line_values = flux_matrix.mean(axis=0)
        std_line_values = flux_matrix.std(axis=0)

        # Recover Direct method + HIIchem fitting
        srPlus_lineLabels = bayes_inputs_dict['Input_data']['lineLabels_list']
        inFluxPlus = bayes_inputs_dict['Input_data']['inputFlux_array']
        inErrPlus = bayes_inputs_dict['Input_data']['inputErr_array']
        flux_matrix_plus = bayes_plus_PI_pickle['trace']['calcFluxes_Op']
        mean_plus_line_values = flux_matrix_plus.mean(axis=0)
        std_plus_line_values = flux_matrix_plus.std(axis=0)

        # Recover Photo-ionization Bayesian fitting
        pI_lineLabels = pI_inputs_dict['Input_data']['lineLabels_list']
        inFlux_PI = pI_inputs_dict['Input_data']['inputFlux_array']
        inErr_PI = pI_inputs_dict['Input_data']['inputErr_array']

        flux_tensor_dict = {}
        trace = pI_pickle['trace']
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

        txt_headers = ['Observed_flux',
                       'Observed_err',
                       'DirectMethod_fitFlux',
                       'DirectMethod_fitSigma',
                       'DirectMethod_difPercentage',
                       'HIICHImistry_Spherical_fitFlux',
                       'HIICHImistry_Spherical_fitSigma',
                       'HIICHImistry_Spherical_difPercentage',
                       'HIICHImistry_PlaneParallel_fitFlux',
                       'HIICHImistry_PlaneParallel_fitSigma',
                       'HIICHImistry_PlaneParallel_difPercentage',
                       'BayesHIICHImistry_fitFlux',
                       'BayesHIICHImistry_fitSigma',
                       'BayesHIICHImistry_difPercentage',
                       'DirectMethod-BayesHIICHImistry_fitFlux',
                       'DirectMethod-BayesHIICHImistry_fitSigma',
                       'DirectMethod-BayesHIICHImistry_difPercentage']
        tableDF = pd.DataFrame(columns=txt_headers)

        pdf = PdfPrinter()
        pdf.create_pdfDoc(pdf_type=None)
        pdf.pdf_insert_table(table_headers, addfinalLine=True)

        for i, lineLabel in enumerate(tableLines):

            row_data = ['-'] * len(table_headers)
            row_txt = ['-'] * (len(txt_headers))

            row_data[0] = latexLabel_array[i]

            obsFlux, obsErr = objLinesDF.loc[lineLabel, 'obsFlux'], objLinesDF.loc[lineLabel, 'obsFluxErr']
            row_data[1] = r'${:.3f}\pm{:.3f}$'.format(obsFlux, obsErr)
            row_txt[0], row_txt[1] = f'{obsFlux:.3f}', f'{obsErr:.3f}'

            # Bayesian direct method fitting
            if lineLabel in sr_lineLabels:
                i_line = sr_lineLabels.index(lineLabel)
                # row_data[1] =  r'${:.3f}\pm{:.3f}$'.format(inFlux[i_line], inErr[i_line])
                row_data[2] = r'${:.3f}\pm{:.3f}$'.format(mean_line_values[i_line], std_line_values[i_line])
                row_data[2] += percent_func(obsFlux, mean_line_values[i_line])

                row_txt[2], row_txt[3] = f'{mean_line_values[i_line]:.3f}', f'{std_line_values[i_line]:.3f}'
                row_txt[4] = percent_func_txt(obsFlux, mean_line_values[i_line])

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

                row_txt[5], row_txt[6] = f'{lineFlux:.3f}', f'{lineFluxErr:.3f}'
                row_txt[7] = percent_func_txt(obsFlux, lineFlux)

            section = f'HII_Tef_fit_C17_bb_Teff_30-90_pp.dat_{cycle}_sed'
            if lineLabel in measurements_dict[section]:
                lineInt, lineIntErr = np.array(measurements_dict[section][lineLabel], dtype=float)
                lineFlux, lineFluxErr = lineInt / np.power(10, cHbeta[0] * flambda), lineIntErr / np.power(10, cHbeta[0] * flambda)
                row_data[4] = r'${:.3f}\pm{:.3f}$'.format(lineFlux, lineFluxErr)
                row_data[4] += percent_func(obsFlux, lineFlux)

                row_txt[8], row_txt[9] = f'{lineFlux:.3f}', f'{lineFluxErr:.3f}'
                row_txt[10] = percent_func_txt(obsFlux, lineFlux)

            if lineLabel in flux_tensor_dict:
                lineFlux, lineFluxErr = np.mean(flux_tensor_dict[lineLabel]), np.std(flux_tensor_dict[lineLabel])
                row_data[5] = r'${:.3f}\pm{:.3f}$'.format(lineFlux, lineFluxErr)
                row_data[5] += percent_func(obsFlux, lineFlux)

                row_txt[11], row_txt[12] = f'{lineFlux:.3f}', f'{lineFluxErr:.3f}'
                row_txt[13] = percent_func_txt(obsFlux, lineFlux)

            if lineLabel in srPlus_lineLabels:
                i_line = sr_lineLabels.index(lineLabel)
                row_data[6] = r'${:.3f}\pm{:.3f}$'.format(mean_plus_line_values[i_line], std_plus_line_values[i_line])
                row_data[6] += percent_func(obsFlux, mean_plus_line_values[i_line])

                row_txt[14], row_txt[15] = f'{mean_plus_line_values[i_line]:.3f}', f'{std_plus_line_values[i_line]:.3f}'
                row_txt[16] = percent_func_txt(obsFlux, mean_plus_line_values[i_line])

            pdf.addTableRow(row_data, last_row=False)
            tableDF.loc[lineLabel] = row_txt

        pdf.table.add_hline()

        table_file = tables_folder/f'{obj}_lineFitFluxComparison'
        pdf.generate_pdf(table_file, clean_tex=True)

        # Save the table as a dataframe.
        with open(f'{table_file}.txt', 'wb') as output_file:
            string_DF = tableDF.to_string()
            output_file.write(string_DF.encode('UTF-8'))

