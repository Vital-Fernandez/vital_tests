import numpy as np
from pathlib import Path
import src.specsiser as sr

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

    if i == 2:

        # Declare input files
        objFolder = resultsFolder / f'{obj}'
        fits_file = dataFolder / f'{obj}_{ext}.fits'
        results_file = objFolder / f'{obj}_{ext}_measurements.txt'
        lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'
        fitResults_file = objFolder / f'{obj}_{ext}_fitting_{cycle}.txt'
        simConf = dataFolder / f'{obj}_config.txt'

        # Declare output files
        outputDb = objFolder / f'{obj}_{ext}_pI_fitting_{cycle}.db'
        outputTxt = objFolder / f'{obj}_{ext}_pI_fitting_{cycle}.txt'

        # Load data
        objParams = sr.loadConfData(simConf, group_variables=False)
        results_dict = sr.loadConfData(results_file, group_variables=False)

        # linesAnalysis = ['O2_3726A_m', 'He1_4471A', 'He1_5876A', 'He2_4686A']
        objLinesDF = sr.import_emission_line_data(lineLog_file, include_lines=objParams[obj]['input_lines'])
        fit_results_dict = sr.loadConfData(fitResults_file, group_variables=False)

        abund_dict = fit_results_dict['Fitting_results']
        O2_abund = np.power(10, abund_dict['O2'] - 12)
        O3_abund = np.power(10, abund_dict['O3'] - 12)

        OH = np.log10(O2_abund[0] + O3_abund[0]) + 12
        OH_err = np.log10(np.sqrt(O2_abund[1] ** 2.0 + O3_abund[1] ** 2.0)) + 12

        cHbeta = fit_results_dict['Fitting_results']['cHbeta']

        # Declare sampler
        obj1_model = sr.SpectraSynthesizer()

        # Declare simulation physical properties
        objRed = sr.ExtinctionModel(Rv=objParams['simulation_properties']['R_v'],
                                    red_curve=objParams['simulation_properties']['reddenig_curve'],
                                    data_folder=objParams['data_location']['external_data_folder'])

        # Declare region physical model
        obj1_model.define_region(objLinesDF, extinction_model=objRed)

        # Replace the flamda values
        idcs_lines = objLinesDF.index.isin(obj1_model.lineLabels)
        flambda_new = objLinesDF.loc[idcs_lines, 'f_lambda']
        obj1_model.lineFlambda = flambda_new

        # Declare sampling properties
        obj1_model.simulation_configuration(objParams['inference_model_configuration']['parameter_list'],
                                            prior_conf_dict=objParams['priors_configuration'],
                                            photo_ionization_grid=True)

        # Declare simulation inference model
        obj1_model.inference_photoionization(OH=OH, cHbeta=cHbeta[0], OH_err=OH_err)

        print(f'- Fitting: {obj}, at OH = {OH}, cHbeta = {cHbeta[0]}')

        # Run the simulation
        obj1_model.run_sampler(objFolder/outputDb, 5000, 2000, njobs=1)

        # Plot the results
        fit_results = sr.load_MC_fitting(outputDb)

        print('-- Model parameters table')
        figure_file = objFolder / f'{obj}_pI_fitting_MeanOutputs'
        obj1_model.table_mean_outputs(figure_file, fit_results)

        print('-- Model parameters posterior diagram')
        figure_file = objFolder / f'{obj}_pI_fitting_ParamsPosteriors.png'
        obj1_model.tracesPosteriorPlot(figure_file, fit_results)

        print('-- Model parameters corner diagram')
        figure_file = objFolder / f'{obj}_pI_fitting_cornerPlot.png'
        obj1_model.corner_plot(figure_file, fit_results)

        print('-- Model emission flux posteriors')
        figure_file = objFolder/f'{obj}_pI_EmFluxPosteriors.png'
        obj1_model.fluxes_photoIonization_distribution(figure_file, fit_results, combined_dict={'O2_3726A_m': 'O2_3726A-O2_3729A',
                                                                                'S2_6716A_m': 'S2_6716A-S2_6731A'})

        print('-- Model emission flux table')
        figure_file = objFolder/f'{obj}_pI_EmFluxPosteriors'
        obj1_model.table_line_fluxes_photoIoniz(figure_file, fit_results, combined_dict={'O2_3726A_m': 'O2_3726A-O2_3729A',
                                                                                'S2_6716A_m': 'S2_6716A-S2_6731A'})