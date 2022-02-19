import src.specsiser as sr
from pathlib import Path
from delete.data_printing import PdfPrinter

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

papaderos_fittings_folder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/Papaderos_Full/out_FD')
ssp_lib_list = ['ChabP94Z5N295', 'midChab94', 'Z3SalpP00', 'Z4P00']
"""
_1D: One dimensional spectra contains:
_EL: Emission-line results
_ST: Statistics: By-products from FADO
_DE: Differential evolution parameters
"""
ext_fado_fits_list = ['1D', 'EL', 'ST', 'DE']
params_per_ext = {'1D':['ARQ_BASE', 'ARQ_CONF', 'CHI2_VAL', 'CHI2_RED', 'NUMPARAM', 'NUM_BASE'],
                  'EL':['NE_LINES', 'N_PARAMS', 'LAMBDA_0', 'GALSNORM', 'TELECTRO', 'DELECTRO', 'GEXTINCT', 'GEXTBDEV', 'GNEBULAR', 'GNEBBDEV'],
                  'ST':['GEXTINCT', 'GNEBULAR', 'AVE_LAGE', 'DEV_LAGE', 'AVE_MAGE', 'DEV_MAGE', 'AVELLAGE', 'DEVLLAGE', 'AVELMAGE', 'DEVLMAGE',
                        'AVE_LMET', 'DEV_LMET', 'AVE_MMET', 'DEV_MMET',
                        'LOGPEAVE', 'LOGPEDEV', 'LOGPCAVE', 'LOGPCDEV',
                        'LOGMCAVE', 'LOGMCDEV', 'LOGPEAVE', 'LOGPEDEV',
                        'LOGMEAVE', 'LOGMEDEV'],
                  'DE':[]}

ext = 'BR'
cycle = 'it3'

dict_ref = {'Library label': 'ARQ_BASE',
            'Number bases': 'NUM_BASE',
            r'$\chi^2$': 'CHI2_RED',
            r'$T_{e} (K)$': 'TELECTRO',
            r'$n_{e} (cm^{-3})$': 'DELECTRO',
            r'$A_{V, neb}$': 'GNEBULAR',
            r'$A_{V, \star}$': 'GEXTINCT',
            'Mass ever formed': 'LOGMEAVE',
            'Mass presently availabled': 'LOGMEAVE',
            'Mass > 1 Myr': 'LOGPEAVE'}

row_headers = list(dict_ref.keys())
tables_folder = Path(obsData['file_information']['tables_folder'])

for i, obj in enumerate(objList):

    pdf = PdfPrinter()
    row_object = [''] * len(row_headers)
    row_object[5] = obj
    pdf.create_pdfDoc(pdf_type='table')
    pdf.pdf_insert_table(row_object, addfinalLine=False)
    pdf.table.add_hline()
    pdf.addTableRow(row_headers, last_row=False, rounddig=2)
    pdf.table.add_hline()

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'

    # Declare output file
    table_file = tables_folder/f'{obj}_FadoFitting'

    # # Load the data
    results_dict = sr.loadConfData(results_file, group_variables=False)

    for j, ssp_lib in enumerate(ssp_lib_list):

        section_label = f'{ssp_lib}_FADO'
        section_results = results_dict[section_label]

        row_data = []
        for column, param_code in dict_ref.items():
            param_value = section_results[param_code]
            row_data.append(param_value)
        row_data[0] = ssp_lib

        pdf.addTableRow(row_data, last_row=False, rounddig=2)

    pdf.table.add_hline()
    pdf.generate_pdf(table_file)
