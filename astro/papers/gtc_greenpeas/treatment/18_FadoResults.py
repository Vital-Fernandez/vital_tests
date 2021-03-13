import src.specsiser as sr
from pathlib import Path

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

for i, obj in enumerate(objList):

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'

    # # Load the data
    # results_dict = sr.loadConfData(results_file, group_variables=False)

    for j, ssp_lib in enumerate(ssp_lib_list):

        # Results container
        fado_measurements = {}
        for k, ext_fit in enumerate(ext_fado_fits_list):

            fits_address = papaderos_fittings_folder/ssp_lib/f'{obj}_FD.cxt.FADO_{ext_fit}.fits'
            data, header = sr.import_fits_data(fits_address, instrument=None, frame_idx=0)

            for param in params_per_ext[ext_fit]:
                print(f'In {ssp_lib} fittig: {param} = {header[param]} ')
                fado_measurements[param] = header[param]

        # Store the results
        print(fado_measurements)
        section_label = f'{ssp_lib}_FADO'
        sr.parseConfDict(results_file, fado_measurements, section_label, clear_section=True)
