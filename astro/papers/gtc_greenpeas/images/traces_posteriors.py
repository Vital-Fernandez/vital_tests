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
        figure_file = objFolder/f'{obj}_DirectMethod_posteriors.png'
        figure_file_flux = objFolder/f'{obj}_DirectMethod_ObsVsFitFluxes.png'
        table_file_fluxes = objFolder/f'{obj}_DirectMethod_ObsVsFitFluxes'
        print(outputDb)
        print(figure_file)

        # Load the results
        fit_pickle = sr.load_MC_fitting(outputDb)
        fit_inputs_dict = sr.loadConfData(outputTxt)
        obj1_model = sr.SpectraSynthesizer()

        # Format outputs from old storing container
        inLines = fit_inputs_dict['Input_data']['lineLabels_list']
        inFluxes, inErr = fit_inputs_dict['Input_data']['inputFlux_array'], fit_inputs_dict['Input_data']['inputErr_array']
        inParameters = np.array(['T_high', 'n_e', 'Ar3', 'Ar4', 'Fe3', 'He1', 'He2', 'N2', 'Ne3', 'O2', 'O3', 'S2', 'S3', 'cHbeta'])
        traces_dict = {param: fit_pickle['trace'][param] for param in inParameters}
        for i, lineLabel in enumerate(inLines):
            traces_dict[lineLabel] = fit_pickle['trace']['calcFluxes_Op'][:, i]

        # Printint results
        # obj1_model.tracesPosteriorPlot(figure_file, inParameters, traces_dict, dark_mode=False)
        # obj1_model.fluxes_distribution(figure_file_flux, inLines, inFluxes, inErr, traces_dict)
        obj1_model.table_line_fluxes(table_file_fluxes, inLines, inFluxes, inErr, traces_dict)


