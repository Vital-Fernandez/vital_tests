import lime
import numpy as np
from astro.papers.muse_CGCG007.muse_CGCG007_methods import save_log_maps, latex_labels, signif_figures, param_units
from lime.plots import PdfMaker
from pathlib import Path

def to_latex_format(array, format_default=True):

    label = r'${}^{{{}}}_{{{}}}$'.format(array[0], array[1], array[2])

    return label


param_list = ['n_e', 'T_low', 'cHbeta', 'Ar4', 'Ar3', 'O2', 'O3', 'N2', 'He1', 'S2', 'S3']
header_list = ['Parameter', 'All voxels', 'Mask 0 (11 voxels)', 'Mask 1 (91 voxels)', 'Mask 2 (382 voxels)']
mask_list = ['MASK_0', 'MASK_1', 'MASK_2']

obsData = lime.load_cfg('../muse_CGCG007.ini')

tablesFolder = Path(obsData['data_location']['tables_folder'])

pdf = PdfMaker()
pdf.create_pdfDoc(pdf_type='table')
pdf.pdf_insert_table(header_list, addfinalLine=True)

for i, param in enumerate(param_list):

    row = ['-'] * len(header_list)

    row[0] = latex_labels[param]

    row[1] = to_latex_format(obsData[f'Global_direct_method'][f'{param}_array'])

    for i_mask, mask_name in enumerate(mask_list):
        row[i_mask + 2] = to_latex_format(obsData[f'{mask_name}_direct_method'][f'{param}_array'])

    pdf.addTableRow(row, last_row=False)

pdf.table.add_hline()
pdf.generate_pdf(tablesFolder/'direct_method_ionic_abundances')

