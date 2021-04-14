import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import compute_arms_flambda, deredd_fluxes, normalize_flux, table_fluxes
import pyneb as pn
from src.specsiser.data_printing import PdfPrinter, latex_labels, numberStringFormat
import pylatex

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)
objList = obsData['file_information']['object_list']
resultsFolder = Path(obsData['file_information']['results_folder'])
tables_folder = Path(obsData['file_information']['tables_folder'])

ext = 'BR'
cycle = 'it3'

sample_Results = {}

# Personal order:
rows = ['n_e', 'T_high', 'cHbeta', 'Ar3', 'Ar4', 'Fe3', 'N2', 'Ne3', 'O2', 'O3', 'S2', 'S3', 'He1', 'He2']

# Load the data from every object
for i, obj in enumerate(objList[:3]):

    objFolder = resultsFolder / f'{obj}'
    outputDb = objFolder / f'{obj}_{ext}_fitting_{cycle}.db'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'

    fit_results = sr.load_MC_fitting(outputDb)
    results_dict = sr.loadConfData(results_file, group_variables=False)

    # Prepare the data
    total_params_list = np.array(list(fit_results['Fitting_results'].keys()))
    traces_param = []
    for param in total_params_list:
        if '_Op' not in param:
            if param in fit_results['trace'].varnames:
                traces_param.append(param)
    traces_param = np.array(traces_param)

    obj_Results = {}
    traces = fit_results['trace']
    for trace_code in traces_param:

        trace_array = traces[trace_code]
        mean_value = np.mean(trace_array)
        std_dev = np.std(trace_array)
        obj_Results[trace_code] = (mean_value, std_dev)

    sample_Results[obj] = obj_Results

# Generate the table
pdf = PdfPrinter()
pdf.create_pdfDoc(pdf_type=None)
pdf.pdf_insert_table(['Parameter', 'SHOC148', 'GP101157', 'GP121903'])

for i, param in enumerate(rows):
    row_data = ['-'] * 4
    row_data[0] = latex_labels[param]

    for j, obj in enumerate(objList[:3]):
        param_mean, param_std = sample_Results[obj][param]

        if param != 'He2':
            round_n = 0 if param_mean > 10 else 2
            param_formated = r'${}\pm{}$'.format(numberStringFormat(param_mean, round_n), numberStringFormat(param_std, round_n))
        else:
            param_formated = r'$\num{{{:.1e}}}\pm\num{{{:.0e}}}$'.format(param_mean, param_std)
        row_data[j+1] = param_formated

    pdf.addTableRow(row_data, last_row=False)
pdf.table.add_hline()

table_file = tables_folder/f'sample_bayesian_fitting'
pdf.generate_pdf(table_file)