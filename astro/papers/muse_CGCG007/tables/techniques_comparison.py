import lime
import numpy as np
import pylatex
from astro.papers.muse_CGCG007.muse_CGCG007_methods import save_log_maps, latex_labels, signif_figures, param_units
from lime.plots import PdfMaker
from pathlib import Path

def to_latex_format(array, format_default=True):

    label = r'${}^{{{}}}_{{{}}}$'.format(array[0], array[1], array[2])

    return label

# Data location
obsData = lime.load_cfg('../muse_CGCG007.ini')
tablesFolder = Path(obsData['data_location']['tables_folder'])

# Table Header structure
param_list = ['OH', 'NO', 'logU']
mask_list = ['MASK_0', 'MASK_1', 'MASK_2']
method_list = ['_direct_method', '_HII-CHI-mistry', '_GridSampling']
header_list = ['',
               pylatex.MultiColumn(3, align='c', data='All voxels'),
               pylatex.MultiColumn(3, align='c', data='Mask 0 (11 voxels)'),
               pylatex.MultiColumn(3, align='c', data='Mask 1 (91 voxels)'),
               pylatex.MultiColumn(3, align='c', data='Mask 2 (382 voxels)')]
subheader_list = ['Parameter',
                  'Direct method', 'HII-CHI-mistry', 'Neural Sampling',
                  'Direct method', 'HII-CHI-mistry', 'Neural Sampling',
                  'Direct method', 'HII-CHI-mistry', 'Neural Sampling',
                  'Direct method', 'HII-CHI-mistry', 'Neural Sampling']
table_header_format = 'lcccccccccccc'


# Pylatex object
pdf = PdfMaker()
pdf.create_pdfDoc(pdf_type='table')
pdf.pdf_insert_table(header_list, table_format=table_header_format, addfinalLine=False)
pdf.addTableRow(subheader_list, last_row=True)

for i, param in enumerate(param_list):

    row = ['-'] * len(subheader_list)

    row[0] = latex_labels[param]

    for j, method_name in enumerate(method_list):

        results_label = f'Global{method_name}'
        param_label = f'{param}_array'
        if param_label in obsData[results_label]:
            row[j + 1] = to_latex_format(obsData[results_label][param_label])

    i_row = 4
    for i_mask, mask_name in enumerate(mask_list):

        for i_method, method_name in enumerate(method_list):

            results_label = f'{mask_name}{method_name}'
            param_label = f'{param}_array'
            if param_label in obsData[results_label]:
                row[i_row] = to_latex_format(obsData[results_label][param_label])
            i_row += 1

    pdf.addTableRow(row, last_row=False)

pdf.table.add_hline()
pdf.generate_pdf(tablesFolder/'methodology_results')

