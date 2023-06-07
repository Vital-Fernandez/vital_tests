import numpy as np
import pandas as pd
import src.specsiser as sr

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import read_lines_fits, compute_cHbeta, safe_cfg
from lime.tools import label_decomposition
from lime import load_cfg, save_line_log

from fastprogress import fastprogress
fastprogress.printing = lambda: True

# Cfg file
cfg_file = Path(r'D:\Pycharm Projects\vital_tests\astro\papers\muse_CGCG007\muse_CGCG007.ini')
obsCfg = load_cfg(cfg_file)

# Number of lines per fit
nights_range = range(1, 4)

# Reddening parameters
red_curve = obsCfg['Extinction']['red_law']
R_v = obsCfg['Extinction']['R_v']

# Results folder
output_folder = Path(r'D:\Dropbox\Astrophysics\Data\CGCG0707_mike')

# Conf file
chem_conf = output_folder/f'MIKE_CGCG007_chemical_fitting.txt'
objParams = load_cfg(chem_conf)

for i_night in nights_range:

    # Declare sampler
    obj1_model = sr.SpectraSynthesizer()

    input_lines = objParams['inference_model_configuration']['input_lines']
    lines_log_address = output_folder/f'MIKE_CGCG007_linelog_night{i_night}.txt'
    objLinesDF = sr.import_emission_line_data(lines_log_address, include_lines=input_lines)
    output_db = output_folder/f'MIKE_CGCG007_db{i_night}'

    lineLabels = objLinesDF.index.values
    lineIons = objLinesDF['ion'].values
    lineFluxes = objLinesDF['intg_flux'].values
    lineErr = objLinesDF['intg_err'].values

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

    # Declare region physical model
    obj1_model.define_region(lineLabels, lineFluxes, lineErr, objIons, objRed, objChem)

    # Declare sampling properties
    obj1_model.simulation_configuration(objParams['inference_model_configuration']['parameter_list'],
                                        prior_conf_dict=objParams['priors_configuration'],
                                        photo_ionization_grid=False)

    # Declare simulation inference model
    obj1_model.inference_model()

    # Run the simulation
    obj1_model.run_sampler(2500, 4000, nchains=3, njobs=1)
    obj1_model.save_fit(output_db)

    # Load the results
    fit_pickle = sr.load_MC_fitting(output_db)
    inLines, inParameters = fit_pickle['inputs']['lines_list'], fit_pickle['inputs']['parameter_list']
    inFlux, inErr = fit_pickle['inputs']['line_fluxes'].astype(float), fit_pickle['inputs']['line_err'].astype(float)
    traces_dict = fit_pickle['outputs']

    # Print the results
    print('-- Model parameters table')
    figure_file = output_folder/f'MIKE_CGCG007_night{i_night}_MeanOutputs'
    obj1_model.table_mean_outputs(figure_file, inParameters, traces_dict)

    print('-- Flux values table')
    figure_file = output_folder/f'MIKE_CGCG007_night{i_night}_FluxComparison'
    obj1_model.table_line_fluxes(figure_file, inLines, inFlux, inErr, traces_dict)

    print('-- Model parameters posterior diagram')
    figure_file = output_folder/f'MIKE_CGCG007_night{i_night}_ParamsPosteriors.png'
    obj1_model.tracesPosteriorPlot(figure_file, inParameters, traces_dict)

    print('-- Line flux posteriors')
    figure_file = output_folder/f'MIKE_CGCG007_night{i_night}_lineFluxPosteriors.png'
    obj1_model.fluxes_distribution(figure_file, inLines, inFlux, inErr, traces_dict)

    print('-- Model parameters corner diagram')
    figure_file = output_folder/f'MIKE_CGCG007_night{i_night}_cornerPlot.png'
    obj1_model.corner_plot(figure_file, inParameters, traces_dict)


'''
Loop through the lineArm dict and put it there (linelabel, flux and err)

Sort it by wavelength

normalized flux by Hbeta

save it as a lines log


'''
