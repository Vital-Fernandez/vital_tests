import src.specsiser as sr
import numpy as np
from pathlib import Path
from src.specsiser.components.starContinuum_functions import SSPsynthesizer
from delete.data_printing import PdfPrinter

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
dataFolder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data')
resultsFolder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/treatment')
starlight_folder = Path('D:/Dropbox/Astrophysics/Tools/Starlight')

obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)
objList = obsData['file_information']['object_list']
objRefNames = obsData['file_information']['refName_list']

ext = '_BR'
cycle = 'it3'
cycle_ref = 'Second_cycle'

tables_folder = Path(obsData['file_information']['tables_folder'])
pdfTableFile = tables_folder/'SSP_results'

row_headers = ['Green Pea Galaxy', '$M_{\star\,total}$', '$M_{\star\,young}$', r'$A_{V}$', r'$c(H\beta)$', r'$S/N$', r'$\chi^{2}_{fit}$']
sub_headers = ['',                  '$log(M_{\odot})$',     r'$\%$',            '',         '',             '',         '']

pdf = PdfPrinter()
pdf.create_pdfDoc(pdf_type='table')
pdf.pdf_insert_table(row_headers, addfinalLine=False)
pdf.addTableRow(sub_headers, last_row=True)

for i, obj in enumerate(objList):

    # Declare files location
    fits_file = dataFolder / f'{obj}{ext}.fits'
    objFolder = resultsFolder / f'{obj}'
    lineLog_file = objFolder / f'{obj}{ext}_linesLog_c2.txt'
    results_file = objFolder / f'{obj}{ext}_measurements.txt'
    objMask = objFolder / f'{obj}{ext}_mask.txt'
    nebCompFile = objFolder / f'{obj}{ext}_NebFlux_{cycle}.txt'
    run_ref = f'{obj}{ext}_{cycle}'
    objSSP_outputFile = starlight_folder/'Output'/f'{obj}{ext}_it2.slOutput'

    # Starlight wrapper
    sw = SSPsynthesizer()

    # Data output
    results_dict = sr.loadConfData(results_file, group_variables=False)
    stellar_fit_dict = results_dict[f'Starlight_run_it2']
    stellar_Wave, obj_input_flux, stellar_flux, fit_output = sw.load_starlight_output(objSSP_outputFile)

    if stellar_fit_dict['A_V_stellarr'] == 0:
        stellar_fit_dict['A_V_stellarr'] = 0.0

    row_data = [objRefNames[i],
                np.round(stellar_fit_dict['Galaxy_mass_Current'], 2),
                np.round(stellar_fit_dict['Galaxy_mass_Percentbelow20Myr'], 2),
                np.round(stellar_fit_dict['A_V_stellarr'], 2),
                np.round(stellar_fit_dict['cHbeta_stellar'], 2),
                np.round(stellar_fit_dict['SN'], 2),
                np.round(stellar_fit_dict['Chi2'], 2)]

    pdf.addTableRow(row_data, last_row=False, rounddig=2)

pdf.table.add_hline()
pdf.generate_pdf(pdfTableFile)
