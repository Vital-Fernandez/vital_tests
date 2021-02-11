from pathlib import Path
import src.specsiser as sr
from src.specsiser.physical_model.chemical_model import Standard_DirectMetchod

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

ext = 'BR'
cycle = 'it3'
MC_steps = 1000

for i, obj in enumerate(objList):

    print(f'- Treating {obj}')

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'
    lineLog_file = objFolder / f'{obj}_{ext}_linesLog_{cycle}.txt'

    # Declare output files
    diagram_file = objFolder / f'{obj}_{ext}_diagChart_{cycle}.png'

    # Load the data
    linesDF = sr.lineslogFile_to_DF(lineLog_file)
    results_dict = sr.loadConfData(results_file, group_variables=False)

    # Extinction parameters
    cHbeta_label = f'cHbeta_{ext}_Hbeta_Hgamma_Hdelta'
    cHbeta = results_dict[f'Extinction_{cycle}'][cHbeta_label]

    # Declare lines to be used in analysis
    idcs_obs = ~linesDF.index.str.contains('_b')
    lineLabels = linesDF.loc[idcs_obs].index
    f_lambda = linesDF.loc[idcs_obs, 'f_lambda'].values
    obsFlux_Err = linesDF.loc[idcs_obs, 'obsFlux':'obsFluxErr'].values
    obsInt_Err = linesDF.loc[idcs_obs, 'obsInt':'obsIntErr'].values

    # Compute matrix fluxes
    cm = Standard_DirectMetchod(n_steps=MC_steps)
    flux_dict = cm.declare_line_fluxes(lineLabels, obsFlux_Err[:, 0], obsFlux_Err[:, 1])

    # Establish extinction correction
    int_dict = cm.red_corr(flux_dict, cHbeta=cHbeta[0], cHbeta_err=cHbeta[1], f_lambda=f_lambda)

    # Plot diagnostics diagram for galaxy
    cm.plot_diag_chart(int_dict, plot_address=diagram_file)

    # Confirm temperature and density diagnostics
    obj_model_conf = obsData[f'{obj}_chemical_model']
    obj_lines_conf = obsData[f'{obj}_line_fitting']
    cm.electron_diagnostics(int_dict,
                            neSII_limit_check=obj_model_conf['nSII_lowerDist_check'],
                            Thigh_diag=obj_model_conf['Te_high_diag'])

    # Compute ionic abundances
    cm.ionic_abundances(int_dict, obj_lines_conf, obj_model_conf)

    # Compute elemental abundances

    # Save results abundances
    cm.save_measurements(results_file, cycle)
