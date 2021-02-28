import numpy as np
from pathlib import Path
import src.specsiser as sr


def check_previous_measurements(objName, parameter_list, measurements_dict, cfg_dict):

    output_dict = {}

    # Chemical abundances
    for param in parameter_list:
        if param in measurements_dict['it3_ionic_Abundances']:
            output_dict[param] = measurements_dict['it3_ionic_Abundances'][param]

    if 'O2_3726A_m' in measurements_dict['it3_ionic_Abundances']:
        output_dict['O2'] = measurements_dict['it3_ionic_Abundances']['O2_3726A_m']

    conversion_dict = {'He1': 'He1r', 'He2': 'He2r'}
    for fit_label, measure_label in conversion_dict.items():
        if measure_label in measurements_dict['it3_ionic_Abundances']:
                output_dict[fit_label] = measurements_dict['it3_ionic_Abundances'][measure_label]

    # Extinction
    cHbeta_key = cfg_dict[objName]['cHbeta_label']
    output_dict['cHbeta'] = measurements_dict['Extinction_it3'][cHbeta_key]

    # Electron parameters
    conversion_dict = {'n_e': 'ne', 'T_low': 'Te_low', 'T_high': 'Te_high'}
    for fit_label, measure_label in conversion_dict.items():
        if measure_label in measurements_dict['it3_electron_parameters']:
            output_dict[fit_label] = measurements_dict['it3_electron_parameters'][measure_label]

    return output_dict


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

    if i == 1:

        # Declare input files
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
        objParams = sr.loadConfData(simConf)
        results_dict = sr.loadConfData(results_file, group_variables=False)
        objLinesDF = sr.import_emission_line_data(lineLog_file, include_lines=objParams['input_lines'])

        # Load previous measurements to compare
        true_values = check_previous_measurements(obj, objParams['parameter_list'], results_dict, obsData)

        # Declare extinction properties
        objRed = sr.ExtinctionModel(Rv=objParams['R_v'],
                                    red_curve=objParams['reddenig_curve'],
                                    data_folder=objParams['external_data_folder'])
        # Declare ion properties
        objIons = sr.IonEmissivity(tempGrid=objParams['temp_grid'],
                                   denGrid=objParams['den_grid'])

        # Generate interpolator from the emissivity grids
        ionDict = objIons.get_ions_dict(np.unique(objLinesDF.ion.values))
        objIons.computeEmissivityGrids(objLinesDF, ionDict, combined_dict=combined_line_dict)

        # Declare chemical model
        objChem = sr.DirectMethod(linesDF=objLinesDF, highTempIons=objParams['high_temp_ions_list'])

        # Declare region physical model
        obj1_model = sr.SpectraSynthesizer()
        obj1_model.define_region(objLinesDF, objIons, objRed, objChem)

        # Replace the flamda values
        idcs_lines = objLinesDF.index.isin(obj1_model.lineLabels)
        flambda_new = objLinesDF.loc[idcs_lines, 'f_lambda']
        # for i, lineLabel in enumerate(obj1_model.lineLabels):
        #     quote = f'{lineLabel} ({lineLabel == objLinesDF.index.values[i]}): {obj1_model.lineFlambda[i]} {flambda_new[i]}'
        #     print(quote)

        obj1_model.lineFlambda = flambda_new

        # Declare sampling properties
        obj1_model.simulation_configuration(objParams['parameter_list'], prior_conf_dict=objParams)

        # Declare simulation inference model
        obj1_model.inference_model(include_Thigh_prior=objParams['T_high_check'])

        # Run the simulation
        obj1_model.run_sampler(objFolder/outputDb, 5000, 2000, njobs=1)

        # Plot the results
        fit_results = sr.load_MC_fitting(outputDb)

        print('-- Model parameters table')
        figure_file = objFolder / f'{obj}_MeanOutputs'
        obj1_model.table_mean_outputs(figure_file, fit_results, true_values=true_values)

        print('-- Flux values table')
        figure_file = objFolder / f'{obj}_FluxComparison'
        obj1_model.table_line_fluxes(figure_file, fit_results)

        print('-- Model parameters posterior diagram')
        figure_file = objFolder / f'{obj}_ParamsPosteriors.png'
        obj1_model.tracesPosteriorPlot(figure_file, fit_results, true_values=true_values)

        print('-- Line flux posteriors')
        figure_file = objFolder / f'{obj}_lineFluxPosteriors.png'
        obj1_model.fluxes_distribution(figure_file, fit_results)

        print('-- Model parameters corner diagram')
        figure_file = objFolder / f'{obj}_cornerPlot.png'
        obj1_model.corner_plot(figure_file, fit_results, true_values=true_values)


