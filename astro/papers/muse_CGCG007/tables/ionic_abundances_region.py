import lime
import numpy as np
from astro.papers.muse_CGCG007.muse_CGCG007_methods import save_log_maps, latex_labels, signif_figures, param_units
from lime.plots import PdfMaker
from pathlib import Path

def to_latex_format(array, format_default=True):

    label = r'${}^{{{}}}_{{{}}}$'.format(array[0], array[1], array[2])

    return label


def to_latex_format_with_uncertainty(array, format_default=True):

    label = r'${}^{{{}}}_{{{}}}\,({})$'.format(array[0], array[1], array[2], array[3])

    return label

param_list = ['n_e', 'T_low', 'cHbeta', 'Ar4', 'Ar3', 'O2', 'O3', 'N2', 'He1', 'S2', 'S3']
header_list = ['Parameter', 'All voxels', 'Region 0 (11 voxels)', 'Region 1 (91 voxels)', 'Region 2 (382 voxels)']
mask_list = ['MASK_0', 'MASK_1', 'MASK_2']

obsData = lime.load_cfg('../muse_CGCG007.ini')

tablesFolder = Path(obsData['data_location']['tables_folder'])

pdf = PdfMaker()
pdf.create_pdfDoc(pdf_type=None)
pdf.pdf_insert_table(header_list, addfinalLine=True)

for i, param in enumerate(param_list):

    row = ['-'] * len(header_list)

    row[0] = latex_labels[param]

    param_value = list(obsData[f'Global_direct_method'][f'{param}_array'])
    param_error = list(obsData[f'Global_err_direct_method'][f'{param}_array'])
    entry_vector = [param_value[0], param_value[1], param_value[2], param_error[0]]
    row[1] = to_latex_format_with_uncertainty(entry_vector)

    for i_mask, mask_name in enumerate(mask_list):
        param_value = list(obsData[f'{mask_name}_direct_method'][f'{param}_array'])
        param_error = list(obsData[f'{mask_name}_err_direct_method'][f'{param}_array'])
        entry_vector = [param_value[0], param_value[1], param_value[2], param_error[0]]
        row[i_mask + 2] = to_latex_format_with_uncertainty(entry_vector)

    pdf.addTableRow(row, last_row=False)

pdf.table.add_hline()
pdf.generate_pdf(tablesFolder/'direct_method_ionic_abundances')

