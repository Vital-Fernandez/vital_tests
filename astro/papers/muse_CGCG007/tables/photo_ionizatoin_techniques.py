import lime
import numpy as np
from pylatex import MultiRow, NoEscape
from astro.papers.muse_CGCG007.muse_CGCG007_methods import save_log_maps, latex_labels, signif_figures, param_units
from lime.plots import PdfMaker
from pathlib import Path


def to_latex_format(array, format_default=True):

    label = r'${}^{{{}}}_{{{}}}$'.format(array[0], array[1], array[2])

    return label


def to_latex_format_with_uncertainty(array, format_default=True):

    label = r'${}^{{{}}}_{{{}}}\,({})$'.format(array[0], array[1], array[2], array[3])

    return label


# Data location
obsData = lime.load_cfg('../muse_CGCG007.ini')
tablesFolder = Path(obsData['data_location']['tables_folder'])

# Table structure
header_list = ['Methodology', 'Parameter', 'All voxels', 'Region 0 (11 voxels)', 'Region 1 (91 voxels)', 'Region 2 (382 voxels)']
model_fitting_params = ['OH', 'NO', 'logU']
direct_method_params = ['OH', 'NO', 'NH', 'ArH', 'SH', 'ICF_S4', 'SO']
methods_list = ['Direct Method', 'HII-CHI-mistry', 'GridSampling', 'GridSampling_MaxErr']
table_header_format = 'c' * len(header_list)
entries_list = ['Global', 'MASK_0', 'MASK_1', 'MASK_2']
convert_grid_names = {'GridSampling': r'\makecell{Neural model fitting \\ (Line error)}',
                      'GridSampling_MaxErr': r'\makecell{Neural model fitting \\ (Maximum line error)}'}


# Pylatex object
pdf = PdfMaker()
pdf.create_pdfDoc(pdf_type=None)
pdf.pdf_insert_table(header_list, table_format=table_header_format)
# pdf.addTableRow(subheader_list, last_row=True)

for i, technique in enumerate(methods_list):

    # Direct method
    if technique == 'Direct Method':

        for j, param in enumerate(direct_method_params):

            row = ['-'] * len(header_list)
            row[0] = MultiRow(len(direct_method_params), data=methods_list[i]) if j == 0 else ''
            row[1] = latex_labels[param]

            for z, region in enumerate(entries_list):
                value = obsData[f'Total_abundances_direct_method'][f'{param}_array']

                entry_vector = [value[0], value[1], value[2]]
                row[z + 2] = to_latex_format(entry_vector)

            pdf.addTableRow(row, last_row=False)

        pdf.table.add_hline()
        pdf.table.add_hline()

    # HII-CHI-mistry results
    if technique == r'HII-CHI-mistry':

        for j, param in enumerate(model_fitting_params):

            row = ['-'] * len(header_list)
            row[0] = MultiRow(3, data=NoEscape(r'\textsc{HII-CHI-mistry}')) if j == 0 else ''
            row[1] = latex_labels[param]

            for z, region in enumerate(entries_list):
                value = obsData[f'{region}_HII-CHI-mistry'][f'{param}_array']
                error = obsData[f'{region}_err_HII-CHI-mistry'][f'{param}_array']

                entry_vector = [value[0], value[1], value[2], error[0]]
                row[z + 2] = to_latex_format_with_uncertainty(entry_vector)

            pdf.addTableRow(row, last_row=False)

        pdf.table.add_hline()
        pdf.table.add_hline()

    if technique in ['GridSampling', 'GridSampling_MaxErr']:

        for j, param in enumerate(model_fitting_params):

            row = ['-'] * len(header_list)
            row[0] = MultiRow(3, data=NoEscape(convert_grid_names[technique])) if j == 0 else ''
            row[1] = latex_labels[param]

            for z, region in enumerate(entries_list):
                value = obsData[f'{region}_{technique}'][f'{param}_array']
                error = obsData[f'{region}_err_{technique}'][f'{param}_array']

                entry_vector = [value[0], value[1], value[2], error[0]]
                row[z + 2] = to_latex_format_with_uncertainty(entry_vector)

            pdf.addTableRow(row, last_row=False)

        if i < len(methods_list) - 1:
            pdf.table.add_hline()
            pdf.table.add_hline()

pdf.table.add_hline()
pdf.generate_pdf(tablesFolder/'methodology_results')

