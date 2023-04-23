import numpy as np
import lime
from pathlib import Path
from astropy.io import fits
from src.specsiser.components.extinction_model import ExtinctionModel

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
tablesFolder = Path(obsData['data_location']['tables_folder'])

# Sample properties
z_objs = obsData['sample_data']['z_array']

# Object
i = 0
obj = objList[i]
objFolder = resultsFolder / obj
fitsLog_address = objFolder / f'{obj}_linesLog.fits'

# Line data
voxel_cords = (167, 167)
log = lime.load_lines_log(fitsLog_address, ext=f'{voxel_cords[0]}-{voxel_cords[1]}_LINELOG')

# Extinction model
Te, ne = 10000.0, 100.0
red_model = ExtinctionModel(Rv=obsData['Extinction']['R_v'], red_curve=obsData['Extinction']['red_law'])
ext_lines = np.array(['H1_4861A', 'H1_6563A', 'H1_8502A', 'H1_8545A', 'H1_8598A', 'H1_8665A', 'H1_8750A',
                     'H1_8863A', 'H1_9015A', 'H1_9229A'])
cHbeta, cHbeta_err = red_model.cHbeta_from_log(log, line_labels=ext_lines, plot_address=False)

norm_label = 'H1_4861A'
flux_norm = log.loc[norm_label, 'gauss_flux'] / 1000.0

# Indeces of blended lines
idcs_nonblended = log['blended_label'] == 'None'

# Labels
log.loc[idcs_nonblended, 'label_i'] = log.loc[idcs_nonblended, 'latexLabel']
log.loc[~idcs_nonblended, 'label_i'] = log.loc[~idcs_nonblended, 'latexLabel'].str[:-1] + '_g$'

# Normalized fluxes
log.loc[idcs_nonblended, 'F_i'] = log.loc[idcs_nonblended, 'intg_flux'] / flux_norm
log.loc[~idcs_nonblended, 'F_i'] = log.loc[~idcs_nonblended, 'gauss_flux'] / flux_norm

# Normalized fluxes
log.loc[idcs_nonblended, 'Err_i'] = log.loc[idcs_nonblended, 'intg_err'] / flux_norm
log.loc[~idcs_nonblended, 'Err_i'] = log.loc[~idcs_nonblended, 'gauss_err'] / flux_norm

# Normalized intensities
int_array = red_model.reddening_correction(log.wavelength.values, log['F_i'].values, log['Err_i'].values,
                                           cHbeta=(cHbeta, cHbeta_err),
                                           R_v=obsData['Extinction']['R_v'],
                                           reddening_curve=obsData['Extinction']['red_law'])

log['I_i'] = int_array[0]
log['I_i_err'] = int_array[1]

# Reddening curve
log['f_lambda'] = red_model.gasExtincParams(log.wavelength.values)
pdfTableFile = tablesFolder / f'example_emission_lines'
table_header_format = 'lccc'
headers = ['Line', r'$f_{\lambda}$', r'$F(\lambda)$', r'$I(\lambda)$']

pdf = lime.io.PdfMaker()
pdf.create_pdfDoc(pdf_type=None)
pdf.pdf_insert_table(headers, table_format=table_header_format, addfinalLine=True)

for j, linelabel in enumerate(log.index):

    label = log.loc[linelabel, 'label_i']
    f_lambda = f'{log.loc[linelabel, "f_lambda"]:0.2f}'
    flux_i, err_f_i = log.loc[linelabel, 'F_i'], log.loc[linelabel, 'Err_i']
    int_i, err_I_i = log.loc[linelabel, 'I_i'], log.loc[linelabel, 'I_i_err']

    flux_entry = r'${:0.0f}\,\pm\,{:0.0f}$'.format(flux_i, err_f_i) if err_f_i > 1.0 else r'${:0.0f}$'.format(flux_i)
    int_entry = r'${:0.0f}\,\pm\,{:0.0f}$'.format(int_i, err_I_i) if err_f_i > 1.0 else r'${:0.0f}$'.format(int_i)

    row_data = [label,
                f_lambda,
                flux_entry,
                int_entry]

    lastRow_check = True if linelabel == log.index[-1] else False
    pdf.addTableRow(row_data, last_row=lastRow_check)

# End table
cHbeta_row = [r'$c(H\beta)$', ''] + [''] * 2
eqwHbeta_row = [r'$-W(\beta)(\AA)$', ''] + [''] * 2
FHbeta_row = [r'$F(H\beta) (10^{-17} \cdot erg\,cm^{-2} s^{-1} \AA^{-1})$', ''] + [''] * 1

idx = i + 2
cHbeta_row[idx] = r'${:0.2f}\,\pm\,{:0.2f}$'.format(cHbeta, cHbeta_err)
eqwHbeta_row[idx] = r'${}\,\pm\,{}$'.format(f'{log.loc[norm_label, "eqw"]:0.0f}', f'{log.loc[norm_label, "eqw_err"]:0.0f}')

F_Hbeta, Err_Hbeta = log.loc[norm_label, "gauss_flux"]/1e-17, log.loc[norm_label, "gauss_err"]/1e-17
FHbeta_row[idx] = r'${}\,\pm\,{}$'.format(f'{F_Hbeta:0.0f}', f'{Err_Hbeta:0.0f}')

pdf.addTableRow(cHbeta_row)
pdf.addTableRow(eqwHbeta_row)
pdf.addTableRow(FHbeta_row, last_row=True)

pdf.generate_pdf(pdfTableFile, clean_tex=False)

# pdf.create_pdfDoc(pdf_type='table')
# pdf.pdf_insert_table(table_headers, addfinalLine=True)
#
# parameter_dict = {'Redshift': z_objs[0],
#                   'Programme ID': '102-B'}
#
# for param, value in parameter_dict.items():
#
#     row = ['-'] * len(table_headers)
#
#     row[0] = param
#     row[1] = value
#
#     pdf.addTableRow(row, last_row=False)
#
# pdf.table.add_hline()
# pdf.generate_pdf(tablesFolder/'CGCG007_properties')
