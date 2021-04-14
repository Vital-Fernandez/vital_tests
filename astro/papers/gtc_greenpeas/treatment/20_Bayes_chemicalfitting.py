import numpy as np
from pathlib import Path
import src.specsiser as sr
from astro.papers.gtc_greenpeas.common_methods import check_previous_measurements

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

ext = 'BR'
cycle = 'it3'
combined_line_dict = {'O2_3726A_m': 'O2_3726A-O2_3729A', 'O2_7319A_m': 'O2_7319A-O2_7330A'}

# Analyse the spectrum
for i, obj in enumerate(objList):

    if i < 3:

        # Declare input files
        print(f'- Treating object ({i}): {obj}')
        objFolder = resultsFolder / f'{obj}'
        fits_file = dataFolder / f'{obj}_{ext}.fits'
        objMask = dataFolder / f'{obj}_{ext}_mask.txt'
        results_file = objFolder / f'{obj}_{ext}_measurements.txt'
        lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'

        # Declare output files
        outputDb = objFolder / f'{obj}_{ext}_fitting_{cycle}.db'
        outputTxt = objFolder / f'{obj}_{ext}_fitting_{cycle}.txt'
        simConf = dataFolder / f'{obj}_config.txt'

        # Load data
        objParams = sr.loadConfData(simConf, group_variables=False)
        objLinesDF = sr.import_emission_line_data(lineLog_file, include_lines=objParams[obj]['input_lines'])

        # Plot the results
        objChem = sr.DirectMethod()
        table_file = objFolder/f'{obj}_elementalabundances'
        objChem.abundances_from_db(outputDb, save_results_address=table_file)

        mean_abund_dict = {}
        for element, trace in objChem.element_traces.items():
            mean_abund_dict[element] = np.array([trace.mean(), trace.std()])

        sr.parseConfDict(results_file, mean_abund_dict, f'Elemental_abundances_{cycle}', clear_section=True)