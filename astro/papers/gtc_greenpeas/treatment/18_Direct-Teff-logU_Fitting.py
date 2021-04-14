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

    if i == 0:

        # Declare input files
        print(f'- Treating object ({i}): {obj}')
        objFolder = resultsFolder / f'{obj}'
        fits_file = dataFolder / f'{obj}_{ext}.fits'
        results_file = objFolder / f'{obj}_{ext}_measurements.txt'
        lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
        fitResults_file = objFolder / f'{obj}_{ext}_fitting_{cycle}.txt'
        simConf = dataFolder / f'{obj}_config.txt'

        # Declare output files
        outputDb = objFolder / f'{obj}_{ext}_Direct-Teff-logU_{cycle}.db'
        outputTxt = objFolder / f'{obj}_{ext}_Direct-Teff-logU_{cycle}.txt'

        # Load data
        objParams = sr.loadConfData(simConf, group_variables=False)
        results_dict = sr.loadConfData(results_file, group_variables=False)
        objLinesDF = sr.import_emission_line_data(lineLog_file, include_lines=objParams[obj]['input_lines'])
        fit_results_dict = sr.loadConfData(fitResults_file, group_variables=False)

        # Load previous measurements to compare
        true_values = check_previous_measurements(obj, objParams['inference_model_configuration']['parameter_list'],
                                                  results_dict, obsData)

        # Declare sampler
        obj1_model = sr.SpectraSynthesizer()

        # Declare simulation physical properties
        objRed = sr.ExtinctionModel(Rv=objParams['simulation_properties']['R_v'],
                                    red_curve=objParams['simulation_properties']['reddenig_curve'],
                                    data_folder=objParams['data_location']['external_data_folder'])

        # Declare ion properties
        objIons = sr.IonEmissivity(tempGrid=objParams['simulation_properties']['temp_grid'],
                                   denGrid=objParams['simulation_properties']['den_grid'])

        # Generate interpolator from the emissivity grids
        ionDict = objIons.get_ions_dict(np.unique(objLinesDF.ion.values))
        objIons.computeEmissivityGrids(objLinesDF, ionDict, combined_dict=combined_line_dict)

        # Declare chemical model
        objChem = sr.DirectMethod(linesDF=objLinesDF, highTempIons=objParams['simulation_properties']['high_temp_ions_list'])

        # Declare region physical model
        obj1_model = sr.SpectraSynthesizer()
        obj1_model.define_region(objLinesDF, objIons, objRed, objChem)

        # Replace the flamda values
        idcs_lines = objLinesDF.index.isin(obj1_model.lineLabels)
        flambda_new = objLinesDF.loc[idcs_lines, 'f_lambda']
        obj1_model.lineFlambda = flambda_new

        # Declare sampling properties
        obj1_model.simulation_configuration(objParams['inference_model_configuration']['parameter_list'],
                                            prior_conf_dict=objParams['priors_configuration'],
                                            photo_ionization_grid=True)

        # Declare simulation inference model
        obj1_model.inference_model(fit_T_low=False)

        # Run the simulation
        obj1_model.run_sampler(objFolder/outputDb, 2000, 2000, njobs=1)

        # Plot the results
        fit_results = sr.load_MC_fitting(outputDb)

        print('-- Model parameters table')
        figure_file = objFolder / f'{obj}_Direct-Teff-logU_MeanOutputs'
        obj1_model.table_mean_outputs(figure_file, fit_results)

        print('-- Model parameters posterior diagram')
        figure_file = objFolder / f'{obj}_Direct-Teff-logU_ParamsPosteriors.png'
        obj1_model.tracesPosteriorPlot(figure_file, fit_results)

        print('-- Model parameters corner diagram')
        figure_file = objFolder / f'{obj}_Direct-Teff-logU_cornerPlot.png'
        obj1_model.corner_plot(figure_file, fit_results)

        print('-- Model emission flux posteriors')
        figure_file = objFolder/f'{obj}_Direct-Teff-logU_EmFluxPosteriors.png'
        obj1_model.fluxes_distribution(figure_file, fit_results, combined_dict=combined_line_dict)

        print('-- Model emission flux table')
        figure_file = objFolder/f'{obj}_Direct-Teff-logU_PhotoFluxPosteriors'
        obj1_model.table_line_fluxes_photoIoniz(figure_file, fit_results, combined_dict={'O2_3726A_m': 'O2_3726A-O2_3729A',
                                                                                'S2_6716A_m': 'S2_6716A-S2_6731A'})