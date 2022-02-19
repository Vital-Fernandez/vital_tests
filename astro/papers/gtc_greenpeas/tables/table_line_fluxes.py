import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import compute_arms_flambda, deredd_fluxes, normalize_flux, table_fluxes
import pyneb as pn
from delete.data_printing import PdfPrinter
import pylatex

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']
tables_folder = Path(obsData['file_information']['tables_folder'])

normFlux = obsData['sample_data']['norm_flux']
w_div_array = obsData['sample_data']['w_div']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']
H1 = pn.RecAtom('H', 1)

dict_linesDF = {}
ext = 'BR'
cycle = 'it3'

for i, obj in enumerate(objList[:3]):

    print(f'\n- Treating {i} :{obj}_{ext}.fits')

    # Declare input files
    fits_file = dataFolder / f'{obj}_{ext}.fits'
    objFolder = resultsFolder / f'{obj}'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
    results_file = objFolder/f'{obj}_{ext}_measurements.txt'

    # Declare output files
    tables_prefix = objFolder/f'{obj}{ext}_linesTable_{cycle}'

    # Load data
    linesDF = sr.lineslogFile_to_DF(lineLog_file)
    results_dict = sr.loadConfData(results_file, group_variables=False)

    # Physical parameters
    Te = results_dict[f'Extinction_{cycle}']['Te_low']
    ne = results_dict[f'Extinction_{cycle}']['ne']
    cHbeta_label = obsData[obj]['cHbeta_label']
    cHbeta = results_dict[f'Extinction_{cycle}'][cHbeta_label]

    # Add new columns to dataframe
    for column in ['obsFlux', 'obsFluxErr', 'f_lambda', 'obsInt', 'obsIntErr']:
        linesDF[column] = np.nan

    # Scale factor for lines in the red arm
    Halpha_beta_ratio = H1.getEmissivity(tem=Te, den=ne, wave=6563) / H1.getEmissivity(tem=Te, den=ne, wave=4861)

    # Distinguish each arm lines
    idcs_blue = linesDF.wavelength < w_div_array[i]
    idcs_red = ~idcs_blue

    # Normalize fluxes blue arm
    obsFlux, obsFluxErr = normalize_flux(linesDF.loc[idcs_blue], norm_line='H1_4861A')
    linesDF.loc[idcs_blue, 'obsFlux':'obsFluxErr'] = np.array([obsFlux, obsFluxErr]).T

    # Normalize fluxes red arm
    obsFlux, obsFluxErr = normalize_flux(linesDF.loc[idcs_red], norm_line='H1_6563A', scale_factor=Halpha_beta_ratio)
    linesDF.loc[idcs_red, 'obsFlux':'obsFluxErr'] = np.array([obsFlux, obsFluxErr]).T

    # Compute reddening law for blue arm
    f_Hbeta, f_blue = compute_arms_flambda(linesDF.loc[idcs_blue], red_law, RV, ref_line='H1_4861A')
    linesDF.loc[idcs_blue, 'f_lambda'] = f_blue

    # Compute reddening law for red arm
    f_Halpha, f_red = compute_arms_flambda(linesDF.loc[idcs_red], red_law, RV, ref_line='H1_6563A')
    linesDF.loc[idcs_red, 'f_lambda'] = f_red

    # Compute line intensities
    obsInt, obsIntErr = deredd_fluxes(linesDF.obsFlux, linesDF.obsFluxErr, cHbeta[0], cHbeta[1], linesDF.f_lambda)
    linesDF.loc[:, 'obsInt':'obsIntErr'] = np.array([obsInt, obsIntErr]).T

    # # Save individual object tables
    idcs_obs = ~linesDF.index.str.contains('_b')
    # table_fluxes(linesDF.loc[idcs_obs].index,
    #              linesDF.loc[idcs_obs].f_lambda,
    #              linesDF.loc[idcs_obs].obsFlux,
    #              linesDF.loc[idcs_obs].obsFluxErr,
    #              linesDF.loc[idcs_obs].obsInt,
    #              linesDF.loc[idcs_obs].obsIntErr,
    #              cHbeta[0],
    #              cHbeta[1],
    #              ref_label='H1_4861A',
    #              ref_flux=linesDF.loc['H1_4861A'].intg_flux,
    #              ref_err=linesDF.loc['H1_4861A'].intg_err,
    #              output_address=tables_prefix)
    #
    # # Save the lines log
    # sr.save_lineslog(linesDF, lineLog_file)

    # Add DF to dict
    dict_linesDF[obj] = linesDF.loc[idcs_obs]


# Global tables
tableDF = pd.DataFrame(columns=['wavelength', 'f_lambda', 'latexLabel', 'blended_label'])
objSubList = ['gp030321', 'gp101157', 'gp121903']

# Combine dictionaries into one DF
for i, obj in enumerate(objSubList):

    # Obj DF
    linesDF_i = dict_linesDF[obj]

    # Add new lines to the master df
    for line in linesDF_i.index:
        if line not in tableDF.index:
            tableDF.loc[line] = [linesDF_i.loc[line, 'wavelength'], linesDF_i.loc[line, 'f_lambda'],
                                 linesDF_i.loc[line, 'latexLabel'], linesDF_i.loc[line, 'blended_label']]
tableDF.sort_values('wavelength', inplace=True)

# Make table
scaleTable = 1000
pdfTableFile, txtTableFile = tables_folder / f'sample_emission_lines', tables_folder / f'sample_emission_lines.txt'
table_header_format = 'lccccccc'
row_headers = ['',
               '',
               pylatex.MultiColumn(2, align='c', data=objSubList[0].upper()),
               pylatex.MultiColumn(2, align='c', data=objSubList[1].upper()),
               pylatex.MultiColumn(2, align='c', data=objSubList[2].upper())]
row_subHeaders = ['Line label', r'$f_{\lambda}$',
                  r'$F(\lambda)$', r'$I(\lambda)$',
                  r'$F(\lambda)$', r'$I(\lambda)$',
                  r'$F(\lambda)$', r'$I(\lambda)$']

# Table heading
pdf = PdfPrinter()
pdf.create_pdfDoc(pdf_type='table')
pdf.pdf_insert_table(row_headers, table_format=table_header_format, addfinalLine=False)
# pdf.table.add_hline(3, 4, cmidruleoption='l{10pt}r{2pt}')
# pdf.table.add_hline(5, 6)
# pdf.table.add_hline(7, 8)
pdf.addTableRow(row_subHeaders, last_row=True)
# Table body
for i, linelabel in enumerate(tableDF.index):

    # Line reference
    if (tableDF.loc[linelabel, 'blended_label'] == 'None') or ('_m' in linelabel):
        label_ref = tableDF.loc[linelabel, 'latexLabel']
    else:
        label_ref = tableDF.loc[linelabel, 'latexLabel'][:-1] + '_g$'

    # Line lambda
    lambda_ref = f'{tableDF.loc[linelabel, "f_lambda"]:0.2f}'

    # List for table row
    row_data = [label_ref, lambda_ref] + ['-'] * 6

    for j, obj in enumerate(objSubList):

        linesDF_i = dict_linesDF[obj]
        if linelabel in linesDF_i.index:
            flux_line = linesDF_i.loc[linelabel, 'obsFlux'] * scaleTable
            flux_err = linesDF_i.loc[linelabel, 'obsFluxErr'] * scaleTable
            int_line = linesDF_i.loc[linelabel, 'obsInt'] * scaleTable
            int_err = linesDF_i.loc[linelabel, 'obsIntErr'] * scaleTable

            row_data[2 + 2 * j] = r'${:0.0f}\,\pm\,{:0.0f}$'.format(flux_line, flux_err)
            row_data[2 + 2 * j + 1] = r'${:0.0f}\,\pm\,{:0.0f}$'.format(int_line, int_err)

    lastRow_check = True if linelabel == tableDF.index[-1] else False
    pdf.addTableRow(row_data, last_row=lastRow_check)

# Ending table
cHbeta_row = [r'$c(H\beta)$', ''] + [''] * 3
eqwHbeta_row = [r'$-W(\beta)(\AA)$', ''] + [''] * 3
FHbeta_row = [r'$F(H\beta) (10^{14} \cdot erg\,cm^{-2} s^{-1} \AA^{-1})$', ''] + [''] * 3
for j, obj in enumerate(objSubList):
    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_BR_measurements.txt'

    results_dict = sr.loadConfData(results_file, group_variables=False)
    linesDF_i = dict_linesDF[obj]
    cHbeta_label = obsData[obj]['cHbeta_label']
    cHbeta = results_dict[f'Extinction_{cycle}'][cHbeta_label]

    lineLog_file_1 = objFolder / f'{obj}_{ext}_linesLog_it1.txt'
    linesDF_1 = sr.lineslogFile_to_DF(lineLog_file_1)
    eqw = linesDF_1.loc["H1_4861A", "eqw"]
    eqw_err =linesDF_1.loc["H1_4861A", "eqw_err"]

    cHbeta = r'${:0.2f}\,\pm\,{:0.2f}$'.format(cHbeta[0], cHbeta[1])
    eqw_Hbeta = r'${}\,\pm\,{}$'.format(f'{eqw:0.0f}', f'{eqw_err:0.0f}')
    F_Hbeta = r'${}$'.format(f'{linesDF_i.loc["H1_4861A", "intg_flux"] / normFlux:0.3f}')

    print(linesDF_i.loc["H1_4861A", "intg_flux"], normFlux, linesDF_i.loc["H1_4861A", "intg_flux"]/normFlux)
    # eqw_Hbeta = f'{linesDF_i.loc["H1_4861A", "eqw"]:0.2f}'
    # F_Hbeta = f'{linesDF_i.loc["H1_4861A", "intg_flux"]/normFlux:0.2f}'

    cHbeta_row[2 + j] = pylatex.MultiColumn(2, align='c', data=pylatex.NoEscape(cHbeta))
    eqwHbeta_row[2 + j] = pylatex.MultiColumn(2, align='c', data=pylatex.NoEscape(eqw_Hbeta))
    FHbeta_row[2 + j] = pylatex.MultiColumn(2, align='c', data=pylatex.NoEscape(F_Hbeta))

pdf.addTableRow(cHbeta_row)
pdf.addTableRow(eqwHbeta_row)
pdf.addTableRow(FHbeta_row, last_row=True)

# Save the pdf table
print(f'-- Saving table at: {pdfTableFile}')
try:
    pdf.generate_pdf(pdfTableFile, clean_tex=False)
except:
    print('-- PDF compilation failure')
