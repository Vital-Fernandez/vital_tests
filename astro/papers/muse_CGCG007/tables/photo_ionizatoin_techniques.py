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
model_fitting_params = ['OH', 'NO', 'logU']
direct_method_params = ['OH', 'NO', 'NH', 'ArH', 'SH', 'ICF_S4', 'SO']


entries_list = ['global', 'MASK_0', 'MASK_1', 'MASK_2']

technique_conf_dict = {'neural_fitting': ['direct_method'],
                       'HII-CHI-mistry': ['HIICHImistry-PopStar'],
                        'GridSampling':  ['localErr', 'HIICHImistry', 'noOII', 'maxErr']}

convert_grid_names = {'direct_method':          'Direct method',
                      'HIICHImistry-PopStar':   NoEscape(r'$\makecell{\textsc{HII-CHI-mistry}}$'),
                      'localErr':               NoEscape(r'\makecell{Neural model fitting \\ (Line error)}'),
                      'HIICHImistry':           NoEscape(r'\makecell{Neural model fitting \\ (HII-CHI-mistry lines and error)}'),
                      'noOII':                  NoEscape(r'\makecell{Neural model fitting \\ (No [OII] lines)}'),
                      'maxErr':                 NoEscape(r'\makecell{Neural model fitting \\ (Uniform maximum error flux)}')}

# Pylatex object
pdf = PdfMaker()
header_list = ['Methodology', 'Parameter', 'All voxels',
               'Region 0 (11 voxels)', 'Region 1 (91 voxels)', 'Region 2 (382 voxels)']
table_header_format = 'c' * len(header_list)
pdf.create_pdfDoc(pdf_type='table')
pdf.pdf_insert_table(header_list, table_format=table_header_format)

for i, technique in enumerate(technique_conf_dict.items()):

    methodology, conf_list = technique

    # Loop through the approaches
    for conf in conf_list:
        print(methodology, conf)

        table_params = direct_method_params if conf == 'direct_method' else model_fitting_params

        for j, param in enumerate(table_params):

            row = ['-'] * len(header_list)
            row[0] = MultiRow(len(table_params), data=convert_grid_names[conf]) if j == 0 else ''
            row[1] = latex_labels[param]

            for z, region in enumerate(entries_list):

                value_array = obsData[f'{methodology}_{conf}'][f'{param}_{region}']

                # Confirm the technique measures that parameter
                entry_vector = [value_array[0], value_array[1], value_array[2]]
                row[z + 2] = to_latex_format(entry_vector)

            pdf.addTableRow(row, last_row=False)

        # Last row
        if (methodology != 'GridSampling') or (methodology == 'GridSampling' and conf != 'maxErr'):
            pdf.table.add_hline()
            pdf.table.add_hline()
        else:
            pdf.table.add_hline()

pdf.generate_pdf(tablesFolder/'methodology_results')

