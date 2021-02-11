import numpy as np
from pathlib import Path
import src.specsiser as sr

def check_previous_measurements(parameter_list, measurements_dict):

        data_dict = {}

        for param in parameter_list:

                if param in ('n_e', 'T_low'):
                        dict_section = measurements_dict['Third_cycle_Ionic_Abundances']
                        value = dict_section[param]
                elif param in ('cHbeta'):
                        dict_section = measurements_dict['Initial_values']
                        value = dict_section['cHbeta_BR_Hbeta_Hgamma_Hdelta']
                else:
                        dict_section = measurements_dict['Third_cycle_Ionic_Abundances']
                        if param in ('Cl3', 'He1r', 'N2', 'O3', 'S2'):
                                if param in dict_section:
                                        value = dict_section[param]
                                else:
                                        value = None

        if param is not None:

                data_dict[param] = value

        return

# # Import the observation data
# obsData = sr.loadConfData('../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini', group_variables=False)
# linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
# data_folder = Path(obsData['file_information']['data_folder'])
# file_list = obsData['file_information']['files_list']
# addressList = list(data_folder/file for file in file_list)

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']
tables_folder = Path(obsData['file_information']['tables_folder'])
normFlux = obsData['sample_data']['norm_flux']


ext = '_BR'
cycle = 'c3'
cycle_ref = 'Second_cycle'
combined_line_dict = {'O2_3726A_m': 'O2_3726A-O2_3729A', 'O2_7319A_b': 'O2_7319A-O2_7330A'}
parameter_list = ['n_e','T_low','cHbeta','Ar3','Ar4','Fe3','He1r','He2r','N2','Ne3','O3','S2','S3']


# Analyse the spectrum
for i, obj in enumerate(objList):

        if i ==2:

                # Declare files location
                fits_file = dataFolder / f'{obj}{ext}.fits'
                objFolder = resultsFolder / f'{obj}'
                results_file = objFolder / f'{obj}{ext}_measurements.txt'
                lineLog_file = objFolder / f'{obj}{ext}_linesLog_{cycle}.txt'
                outputDb = objFolder/f'{obj}{ext}_fitting_{cycle}.db'
                outputTxt = objFolder/f'{obj}{ext}_fitting_{cycle}.txt'
                simConf = dataFolder / f'{obj}_config.txt'
                results_file = objFolder / f'{obj}{ext}_measurements.txt'

                # Load data
                objParams = sr.loadConfData(simConf)
                results_dict = sr.loadConfData(results_file, group_variables=False)
                objLinesDF = sr.import_emission_line_data(lineLog_file, include_lines=objParams['input_lines'])



                # # Recover results from direct method:
                # standard_method_results = dict(n_e = results_dict['Third_cycle_Electron_parameters']['ne'],
                #                                T_low = results_dict['Third_cycle_Electron_parameters']['Te_low'],
                #                                cHbeta = results_dict['[Extinction_c3]']['cHbeta_BR_Halpha_Hbeta_Hgamma_Hdelta'],
                #                                Ar3= results_dict['Third_cycle_Ionic_Abundances']['XXX'],
                #                                Ar4= results_dict['Third_cycle_Ionic_Abundances']['XXX'],
                #                                Fe3= results_dict['Third_cycle_Ionic_Abundances']['XXX'],
                #                                He1r= results_dict['Third_cycle_Ionic_Abundances']['XXX'],
                #                                He2r= results_dict['Third_cycle_Ionic_Abundances']['XXX'],
                #                                N2= results_dict['Third_cycle_Ionic_Abundances']['XXX'],
                #                                Ne3= results_dict['Third_cycle_Ionic_Abundances']['XXX'],
                #                                O3= results_dict['Third_cycle_Ionic_Abundances']['XXX'],
                #                                S2= results_dict['Third_cycle_Ionic_Abundances']['XXX'],
                #                                S3= results_dict['Third_cycle_Ionic_Abundances']['XXX'],)

                # Declare extinction properties
                objRed = sr.ExtinctionModel(Rv=objParams['R_v'],
                                            red_curve=objParams['reddenig_curve'],
                                            data_folder=objParams['external_data_folder'])
                # Declare ion properties
                objIons = sr.IonEmissivity(tempGrid=objParams['temp_grid'],
                                           denGrid=objParams['den_grid'])

                # Generate interpolator from the emissivity grids
                ionDict = objIons.get_ions_dict(np.unique(objLinesDF.ion.values))
                objIons.computeEmissivityGrids(objLinesDF, ionDict, linesDb=sr._linesDb, combined_dict=combined_line_dict)

                # Declare chemical model
                objChem = sr.DirectMethod(linesDF=objLinesDF, highTempIons=objParams['high_temp_ions_list'])

                # Declare region physical model
                obj1_model = sr.SpectraSynthesizer()

                obj1_model.define_region(objLinesDF, objIons, objRed, objChem)

                # Declare sampling properties
                obj1_model.simulation_configuration(objParams['parameter_list'], prior_conf_dict=objParams)

                # Declare simulation inference model
                obj1_model.inference_model(include_Thigh_prior=objParams['T_high_check'])

                # Run the simulation
                obj1_model.run_sampler(objFolder/outputDb, 5000, 2000, njobs=1)

                # Plot the results
                fit_results = sr.load_MC_fitting(outputDb)

                print('-- Model parameters table')
                figure_file = objFolder/f'{obj}_MeanOutputs'
                obj1_model.table_mean_outputs(figure_file, fit_results)

                print('-- Flux values table')
                figure_file = objFolder/f'{obj}_FluxComparison'
                obj1_model.table_line_fluxes(figure_file, fit_results)

                print('-- Model parameters posterior diagram')
                figure_file = objFolder/f'{obj}_ParamsPosteriors.png'
                obj1_model.tracesPosteriorPlot(figure_file, fit_results)

                print('-- Line flux posteriors')
                figure_file = objFolder/f'{obj}_lineFluxPosteriors.png'
                obj1_model.fluxes_distribution(figure_file, fit_results)

                # print('-- Model parameters corner diagram')
                # figure_file = simFolder/f'{objName}_corner'
                # obj1_model.corner_plot(objParams['parameter_list'], traces_dict)
                # obj1_model.savefig(figure_file, resolution=200)



