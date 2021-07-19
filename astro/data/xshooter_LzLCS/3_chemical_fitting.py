from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import src.specsiser as sr

c_KMpS = 299792.458

obsData = sr.loadConfData('./xshooter_LzLCS.ini')
data_folder = Path(obsData['data_location']['data_folder'])
results_folder = Path(obsData['data_location']['results_folder'])
objfile_list = obsData['data_location']['objfile_list']
sigmafile_list = obsData['data_location']['sigmafile_list']
objRef_list = obsData['data_location']['ref_list']
maskfile = obsData['data_location']['generalMask']

wmin_array = obsData['sample_data']['w_min_array']
wmax_array = obsData['sample_data']['w_max_array']
norm_flux = obsData['sample_data']['norm_flux']
z_obj = obsData['sample_data']['z_obj']
profile_conf = obsData['line_fitting']

DF_list = [None, None]
for i, objName in enumerate(objRef_list):

    # input data
    lineslog_file = results_folder/f'{objName}_linesLog.txt'

    # Load data
    linesDF = sr.lineslogFile_to_DF(lineslog_file)
    DF_list[i] = linesDF

# Join the dataframes
objDF = pd.concat(DF_list)

# Chemical configuration
chemical_fit_file = data_folder/'j131037_chemical_conf.txt'
objParams = sr.loadConfData(chemical_fit_file, group_variables=False)

for component in ['narrow', 'wide']:

    # Declare input lines
    input_lines = obsData['chemical_analysis'][f'input_{component}_line_list']
    norm_line = obsData['chemical_analysis'][f'{component}_norm']

    norm_flux = objDF.loc[norm_line, 'gauss_flux']
    idcs_lines = objDF.index.isin(input_lines)
    lineLabels = objDF.loc[idcs_lines].index
    lineIons = objDF.loc[idcs_lines, 'ion'].values
    lineFluxes = objDF.loc[idcs_lines, 'gauss_flux'].values / norm_flux
    lineErr = objDF.loc[idcs_lines, 'gauss_err'].values / norm_flux

    # Declare simulation physical properties
    objRed = sr.ExtinctionModel(Rv=objParams['simulation_properties']['R_v'],
                                red_curve=objParams['simulation_properties']['reddenig_curve'],
                                data_folder=objParams['data_location']['external_data_folder'])

    objIons = sr.IonEmissivity(tempGrid=objParams['simulation_properties']['temp_grid'],
                               denGrid=objParams['simulation_properties']['den_grid'])

    # Generate interpolator from the emissivity grids
    ionDict = objIons.get_ions_dict(np.unique(lineIons))
    objIons.computeEmissivityGrids(lineLabels, ionDict)

    # Declare chemical model
    objChem = sr.DirectMethod(lineLabels, highTempIons=objParams['simulation_properties']['high_temp_ions_list'])

    # Declare sampler
    obj1_model = sr.SpectraSynthesizer()
    obj1_model.define_region(lineLabels, lineFluxes, lineErr, objIons, objRed, objChem)

    # Declare sampling properties
    obj1_model.simulation_configuration(objParams['inference_model_configuration']['parameter_list'],
                                        prior_conf_dict=objParams['priors_configuration'],
                                        photo_ionization_grid=False,
                                        n_regions=1)

    # Declare simulation inference model
    output_db = results_folder/f'J131037_{component}_db'
    obj1_model.inference_model()

    # Run the simulation
    obj1_model.run_sampler(2000, 2000, nchains=4, njobs=4)
    obj1_model.save_fit(output_db)

    # Load the results
    fit_pickle = sr.load_MC_fitting(output_db)
    inLines, inParameters = fit_pickle['inputs']['line_list'], fit_pickle['inputs']['parameter_list']
    inFlux, inErr = fit_pickle['inputs']['line_fluxes'].astype(float), fit_pickle['inputs']['line_err'].astype(float)
    traces_dict = fit_pickle['outputs']

    # Print the results
    #TODO make plots independent of obj1_model
    print('-- Model parameters table')
    figure_file = results_folder/f'{component}_MeanOutputs'
    obj1_model.table_mean_outputs(figure_file, inParameters, traces_dict)

    print('-- Flux values table')
    figure_file = results_folder/f'{component}_FluxComparison'
    obj1_model.table_line_fluxes(figure_file, inLines, inFlux, inErr, traces_dict)

    print('-- Model parameters posterior diagram')
    figure_file = results_folder/f'{component}_ParamsPosteriors.png'
    obj1_model.tracesPosteriorPlot(figure_file, inParameters, traces_dict)

    print('-- Line flux posteriors')
    figure_file = results_folder/f'{component}_lineFluxPosteriors.png'
    obj1_model.fluxes_distribution(figure_file, inLines, inFlux, inErr, traces_dict)