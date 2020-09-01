import numpy as np
from pathlib import Path
import src.specsiser as sr

# Import the observation data
obsData = sr.loadConfData('../gtc_greenpeas_data.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
data_folder = Path(obsData['file_information']['data_folder'])
file_list = obsData['file_information']['files_list']
addressList = list(data_folder/file for file in file_list)

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Establish files location
    objName = obsData['file_information']['object_list'][i]
    fitsFolder, fitsFile = file_address.parent, file_address.name
    lineLogFolder, lineLogFile = fitsFolder / 'flux_analysis', fitsFile.replace('.fits', '_linesLog.txt')
    simFolder, simConf = fitsFolder / 'chemical_analysis',  fitsFile.replace('_BR.fits', '_config.txt')
    inputLinesLog = f'{objName}_inputLinesLog.txt'
    outputDb = f'{objName}_fitting.db'
    outputTxt = f'{objName}_fitting.txt'
    print(f'-{objName}')

    # if objName == 'gp101157': #'gp121903':# 'gp030321'

    obj1_model = sr.SpectraSynthesizer()
    lm = sr.LineMesurerGUI(linesDF_address=lineLogFolder / lineLogFile)

    blended_dict = obsData['blended_groups']
    blended_list = []
    for group in blended_dict:
        blended_list += blended_dict[group].split('-')

    idcs_blended = lm.linesDF.index.isin(blended_list)

    # Asociate the corresponding flux to the appropiate line
    lm.linesDF.insert(loc=1, column='obsFlux', value=np.nan)
    lm.linesDF.insert(loc=2, column='obsFluxErr', value=np.nan)
    flux_Hbeta = lm.linesDF.loc['H1_4861A', 'intg_flux']
    lm.linesDF.loc[idcs_blended, 'obsFlux'] = lm.linesDF.loc[idcs_blended, 'gauss_flux']/flux_Hbeta
    lm.linesDF.loc[idcs_blended, 'obsFluxErr'] = lm.linesDF.loc[idcs_blended, 'gauss_err']/flux_Hbeta
    lm.linesDF.loc[~idcs_blended, 'obsFlux'] = lm.linesDF.loc[~idcs_blended, 'intg_flux']/flux_Hbeta
    lm.linesDF.loc[~idcs_blended, 'obsFluxErr'] = lm.linesDF.loc[~idcs_blended, 'intg_err']/flux_Hbeta
    lm.linesDF.rename(columns={'wavelength': 'obsWave'}, inplace=True)
    lm.save_lineslog(lm.linesDF, simFolder/inputLinesLog)

    # Load simulation parameters
    objParams = sr.loadConfData(simFolder/simConf)

    # Load emission lines
    objLinesDF = sr.import_emission_line_data(simFolder/inputLinesLog, include_lines=objParams['input_lines'])

    # Declare simulation physical properties
    objRed = sr.ExtinctionModel(Rv=objParams['R_v'],
                                red_curve=objParams['reddenig_curve'],
                                data_folder=objParams['external_data_folder'])

    objIons = sr.IonEmissivity(tempGrid=objParams['temp_grid'],
                               denGrid=objParams['den_grid'])

    # Generate interpolator from the emissivity grids
    ionDict = objIons.get_ions_dict(np.unique(objLinesDF.ion.values))
    objIons.computeEmissivityGrids(objLinesDF, ionDict, linesDb=sr._linesDb, combined_dict={'O2_3726A': 'O2_3726A-O2_3729A'})

    # Declare chemical model
    objChem = sr.DirectMethod(linesDF=objLinesDF, highTempIons=objParams['high_temp_ions_list'])

    # Declare region physical model
    obj1_model.define_region(objLinesDF, objIons, objRed, objChem)

    # Declare sampling properties
    obj1_model.simulation_configuration(objParams['parameter_list'], prior_conf_dict=objParams)

    # # Declare simulation inference model
    # obj1_model.inference_emisGrid_model(include_Thigh_prior=objParams['T_high_check'])
    #
    # # Run the simulation
    # obj1_model.run_sampler(simFolder/outputDb, 5000, 2000, njobs=1)

    # Load the results
    obj1_model.load_sampler_results(simFolder/outputDb, simFolder/outputTxt)

    # Plot the results
    traces_dict = sr.load_MC_fitting(simFolder/outputDb, normConstants=objParams)

    # Table mean values
    true_values = None

    print('-- Model parameters table')
    figure_file = simFolder/f'{objName}_MeanOutputs'
    obj1_model.table_mean_outputs(figure_file, traces_dict, true_values)

    print('-- Flux values table')
    figure_file = simFolder/f'{objName}_EmissionFluxComparison'
    obj1_model.table_line_fluxes(figure_file, traces_dict, obj1_model.lineLabels, obj1_model.emissionFluxes,
                                 obj1_model.emissionErr)

    print('-- Model parameters posterior diagram')
    figure_file = simFolder/f'{objName}_ParamsPosteriors'
    obj1_model.tracesPosteriorPlot(objParams['parameter_list'], traces_dict)
    obj1_model.savefig(figure_file, resolution=200)
    #
    # # print('-- Model parameters corner diagram')
    # # figure_file = simFolder/f'{objName}_corner'
    # # obj1_model.corner_plot(objParams['parameter_list'], traces_dict)
    # # obj1_model.savefig(figure_file, resolution=200)
    #
    # print('-- Line flux posteriors')
    # figure_file = simFolder/f'{objName}_lineFluxPosteriors'
    # obj1_model.fluxes_distribution(obj1_model.lineLabels, obj1_model.lineIons, 'calcFluxes_Op', traces_dict,
    #                                obj1_model.emissionFluxes, obj1_model.emissionErr, objLinesDF.latexLabel.values)
    # obj1_model.savefig(figure_file, resolution=200)

