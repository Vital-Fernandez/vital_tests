import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import STANDARD_AXES, DARK_PLOT, background_color, foreground_color
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from src.specsiser.data_printing import latex_labels, numberStringFormat, label_decomposition, PdfPrinter
import matplotlib.gridspec as gridspec
from pylatex import NoEscape


def gen_colorList(vmin=0.0, vmax=1.0, color_palette=None):
    colorNorm = colors.Normalize(vmin, vmax)
    cmap = cm.get_cmap(name=color_palette)
    # return certain color
    # self.cmap(self.colorNorm(idx))

    return colorNorm, cmap


def tracesPosteriorPlot(plot_address, db_dict, true_values=None, plot_conf={}):
    # Prepare the data
    total_params_list = np.array(list(db_dict['Fitting_results'].keys()))
    traces_param = []
    for param in total_params_list:
        if '_Op' not in param:
            traces_param.append(param)
    traces_param = np.array(traces_param)

    # Establish the params analysed
    idx_analysis_lines = np.in1d(traces_param, db_dict['trace'].varnames)
    params_list = traces_param[idx_analysis_lines]

    if true_values is not None:
        trace_true_dict = {}
        for param in params_list:
            if param in true_values:
                trace_true_dict[param] = true_values[param]
    n_traces = len(params_list)

    # Plot format
    size_dict = {'axes.titlesize': 14, 'axes.labelsize': 16, 'legend.fontsize': 12,
                 'xtick.labelsize': 10, 'ytick.labelsize': 10}
    size_dict.update(plot_conf)
    rcParams.update(size_dict)

    fig = plt.figure(figsize=(8, n_traces))
    colorNorm, cmap = gen_colorList(0, n_traces)
    gs = gridspec.GridSpec(n_traces * 2, 4)
    gs.update(wspace=0.2, hspace=1.8)

    traces = db_dict['trace']
    for i in range(n_traces):

        trace_code = params_list[i]
        trace_array = traces[trace_code]

        mean_value = np.mean(trace_array)
        std_dev = np.std(trace_array)

        axTrace = fig.add_subplot(gs[2 * i:2 * (1 + i), :3])
        axPoterior = fig.add_subplot(gs[2 * i:2 * (1 + i), 3])

        # Label for the plot
        if mean_value > 10:
            label = r'{} = ${:.0f}$$\pm${:.0f}'.format(latex_labels[trace_code], mean_value, std_dev)
        else:
            label = r'{} = ${:.3f}$$\pm${:.3f}'.format(latex_labels[trace_code], mean_value, std_dev)

        # # Label for the plot
        # if mean_value > 0.001:
        #     label = r'{} = ${}$ $\pm${}'.format(latex_labels[trace_code], np.round(mean_value, 4), np.round(std_dev, 4))
        # else:
        #     label = r'{} = ${:.3e}$ $\pm$ {:.3e}'.format(latex_labels[trace_code], mean_value, std_dev)

        # Plot the traces
        axTrace.plot(trace_array, label=label, color=cmap(colorNorm(i)))
        axTrace.axhline(y=mean_value, color=cmap(colorNorm(i)), linestyle='--')
        axTrace.set_ylabel(latex_labels[trace_code])

        # Plot the histograms
        axPoterior.hist(trace_array, bins=50, color=cmap(colorNorm(i)), align='left')


        # Plot the axis as percentile
        median, percentile16th, percentile84th = np.median(trace_array), np.percentile(trace_array, 16), np.percentile(
            trace_array, 84)

        # Add true value if available
        if true_values is not None:
            if trace_code in trace_true_dict:
                value_param = trace_true_dict[trace_code]

                # Nominal value and uncertainty
                if isinstance(value_param, (list, tuple, np.ndarray)):
                    nominal_value, std_value = value_param[0], 0.0 if len(value_param) == 1 else value_param[1]
                    axPoterior.axvline(x=nominal_value, color=cmap(colorNorm(i)), linestyle='solid')
                    axPoterior.axvspan(nominal_value - std_value, nominal_value + std_value, alpha=0.5,
                                       color=cmap(colorNorm(i)))

                # Nominal value only
                else:
                    nominal_value = value_param
                    axPoterior.axvline(x=nominal_value, color=cmap(colorNorm(i)), linestyle='solid')

        # Add legend
        axTrace.legend(loc=7)

        # Remove ticks and labels
        if i < n_traces - 1:
            axTrace.get_xaxis().set_visible(False)
            axTrace.set_xticks([])

        axPoterior.yaxis.set_major_formatter(plt.NullFormatter())
        axPoterior.set_yticks([])

        axPoterior.set_xticks([percentile16th, median, percentile84th])
        round_n = 0 if median > 10 else 3
        axPoterior.set_xticklabels(['', numberStringFormat(median, round_n), ''])

        axTrace.set_yticks((percentile16th, median, percentile84th))
        round_n = 0 if median > 10 else 3
        axTrace.set_yticklabels(
            (numberStringFormat(percentile16th, round_n), '', numberStringFormat(percentile84th, round_n)))

    plt.show()

    # if plot_address is not None:
    #     plt.savefig(plot_address, dpi=200, bbox_inches='tight')
    #     plt.close(fig)
    # else:
    #     # plt.tight_layout()
    #     plt.show()

    return


def table_line_fluxes(table_address, db_dict, combined_dict={}, file_type='table'):

    # Table headers
    headers = ['Line', 'Observed flux', 'Mean', 'Standard deviation', 'Median', r'$16^{th}$ $percentil$',
               r'$84^{th}$ $percentil$', r'$Difference\,\%$']

    # Create containers
    tableDF = pd.DataFrame(columns=headers[1:])

    pdf = PdfPrinter()
    pdf.create_pdfDoc(pdf_type=file_type)
    pdf.pdfDoc.append(NoEscape('\definecolor{background}{rgb}{0.169, 0.169, 0.169}'))
    pdf.pdfDoc.append(NoEscape('\definecolor{foreground}{rgb}{0.702, 0.780, 0.847}'))
    pdf.pdfDoc.append(NoEscape(r'\arrayrulecolor{foreground}'))

    pdf.pdf_insert_table(headers, color_font='foreground', color_background='background')

    # Input data
    inputLabels = db_dict['Input_data']['lineLabels_list']
    inFlux, inErr = db_dict['Input_data']['inputFlux_array'], db_dict['Input_data']['inputErr_array']

    # Output data
    flux_matrix = db_dict['trace']['calcFluxes_Op']
    mean_line_values = flux_matrix.mean(axis=0)
    std_line_values = flux_matrix.std(axis=0)
    median_line_values = np.median(flux_matrix, axis=0)
    p16th_line_values = np.percentile(flux_matrix, 16, axis=0)
    p84th_line_values = np.percentile(flux_matrix, 84, axis=0)

    # Array wih true error values for flux
    diff_Percentage = np.round((1 - (median_line_values / inFlux)) * 100, 2)
    diff_Percentage = list(map(str, diff_Percentage))

    ion_array, wave_array, latexLabel_array = label_decomposition(inputLabels, combined_dict=combined_dict)

    for i in range(inFlux.size):

        # label = label_formatting(inputLabels[i])

        row_i = [latexLabel_array[i], inFlux[i], mean_line_values[i], std_line_values[i], median_line_values[i], p16th_line_values[i],
                 p84th_line_values[i], diff_Percentage[i]]

        pdf.addTableRow(row_i, last_row=False if inputLabels[-1] != inputLabels[i] else True, color_font='foreground', color_background='background')
        tableDF.loc[inputLabels[i]] = row_i[1:]

    pdf.generate_pdf(table_address)

    # Save the table as a dataframe.
    with open(f'{table_address}.txt', 'wb') as output_file:
        string_DF = tableDF.to_string()
        output_file.write(string_DF.encode('UTF-8'))


def fluxes_distribution(plot_address, db_dict, n_columns=8, combined_dict={}, plot_conf={}):


    # Input data
    inputLabels = db_dict['Input_data']['lineLabels_list']
    inFlux, inErr = db_dict['Input_data']['inputFlux_array'], db_dict['Input_data']['inputErr_array']
    ion_array, wave_array, latexLabel_array = label_decomposition(inputLabels, combined_dict=combined_dict)
    trace = db_dict['trace']
    model = db_dict['model']

    flux_matrix = db_dict['trace']['calcFluxes_Op']
    median_values = np.median(flux_matrix, axis=0)

    # Declare plot grid size
    n_lines = len(inputLabels)
    n_rows = int(np.ceil(float(n_lines) / float(n_columns)))
    n_cells = n_rows * n_columns

    # Declare figure format
    size_dict = {'figure.figsize': (22, 9), 'axes.titlesize': 14, 'axes.labelsize': 10, 'legend.fontsize': 10,
                 'xtick.labelsize': 8, 'ytick.labelsize': 3}
    size_dict.update(plot_conf)
    rcParams.update(size_dict)

    # self.FigConf(plotSize=size_dict, Figtype='Grid', n_columns=n_columns, n_rows=n_rows)
    fig, axes = plt.subplots(n_rows, n_columns)
    axes = axes.ravel()

    # Generate the color dict
    obsIons = np.unique(ion_array)
    colorNorm, cmap = gen_colorList(0, obsIons.size)
    colorDict = dict(zip(obsIons, np.arange(obsIons.size)))

    # Plot individual traces
    for i in range(n_cells):

        if i < n_lines:

            # Current line
            label = inputLabels[i]
            ion = ion_array[i]
            trace = flux_matrix[:, i]
            median_flux = median_values[i]



            label_mean = 'Mean value: {}'.format(np.around(median_flux, 4))

            norm_trace = (trace-inFlux[i])/trace
            axes[i].hist(norm_trace, histtype='stepfilled', bins=35, alpha=.8, color=cmap(colorNorm(colorDict[ion])),
                         density=False)
            # axes[i].hist(trace, histtype='stepfilled', bins=35, alpha=.8, color=cmap(colorNorm(colorDict[ion])),
            #              density=False)

            label_true = 'True value: {}'.format(np.around(inFlux[i], 3))
            axes[i].axvspan(-inErr[i]/inFlux[i], inErr[i]/inFlux[i], alpha=0.3, color=foreground_color)

            # axes[i].axvline(x=inFlux[i], label=label_true, alpha=0.3, color=foreground_color, linestyle='solid')
            # # axes[i].axvspan(inFlux[i] - inErr[i], inFlux[i] + inErr[i], alpha=0.3, color=foreground_color)
            # axes[i].axvspan(0 - inErr[i], 0 + inErr[i], alpha=0.3, color=foreground_color)

            axes[i].get_yaxis().set_visible(False)
            axes[i].set_yticks([])

            # Plot wording
            axes[i].set_title(latexLabel_array[i])

        else:
            fig.delaxes(axes[i])

    # if plot_address is not None:
    #     plt.savefig(plot_address, dpi=200, bbox_inches='tight')
    #     plt.close(fig)
    # else:
    #     # plt.tight_layout()
    plt.show()

    return


conf_file = Path('/home/vital/PycharmProjects/vital_tests/astro/data/muse/muse_greenpeas.ini')

obsData = sr.loadConfData(conf_file, group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']


folder = Path('/home/vital/Dropbox/Astrophysics/Seminars/UniVapo 2021/')

# Plot set up
defaultConf = DARK_PLOT.copy()
rcParams.update(defaultConf)

# Data location
i = 0
obj = objList[i]
cube_address = fitsFolder / fileList[i]
objFolder = resultsFolder / obj
db_addresss = objFolder / f'{obj}_database.fits'
mask_address = dataFolder / obj / f'{obj}_mask.txt'
folder_voxels = Path('/home/vital/Astro-data/Observations/voxel_data')

idx_j, idx_i = 168, 169

obj1_model = sr.SpectraSynthesizer()

# Plot the results
outputDb = folder_voxels/f'{idx_j}-{idx_i}_fitting.db'
fit_results = sr.load_MC_fitting(outputDb)

print('-- Printing results')
# figure_file = folder_voxels / f'{idx_j}-{idx_i}_ParamsPosteriors.png'
# tracesPosteriorPlot(figure_file, fit_results)

print('-- Printing table')
# figure_file = folder_voxels / f'{idx_j}-{idx_i}_FluxComparison'
# table_line_fluxes(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})

print('-- Printing table')
figure_file = folder_voxels / f'{idx_j}-{idx_i}_lineFluxPosteriors.png'
fluxes_distribution(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})







# figure_file = folder_voxels / f'{idx_j}-{idx_i}_MeanOutputs'
# obj1_model.table_mean_outputs(figure_file, fit_results)

# figure_file = folder_voxels / f'{idx_j}-{idx_i}_FluxComparison'
# obj1_model.table_line_fluxes(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})

# figure_file = folder_voxels / f'{idx_j}-{idx_i}_ParamsPosteriors.png'
# obj1_model.tracesPosteriorPlot(figure_file, fit_results)
#
# figure_file = folder_voxels / f'{idx_j}-{idx_i}_lineFluxPosteriors.png'
# obj1_model.fluxes_distribution(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})
#
# figure_file = folder_voxels / f'{idx_j}-{idx_i}_cornerPlot.png'
# obj1_model.corner_plot(figure_file, fit_results)
