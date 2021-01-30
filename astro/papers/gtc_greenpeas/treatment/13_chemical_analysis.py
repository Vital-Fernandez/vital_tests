import src.specsiser as sr
from pathlib import Path
from src.specsiser.physical_model.chemical_model import Standard_DirectMetchod
from astro.papers.gtc_greenpeas.common_methods import compute_arms_flambda, deredd_fluxes, normalize_flux, table_fluxes
import pyneb as pn

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

ext = '_BR'
cycle = 'c2'
cycle_ref = 'Second_cycle'

counter = 0
for i, obj in enumerate(objList):

    if i > 2:

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objFolder = resultsFolder/f'{obj}'
        lineLog_file = objFolder/f'{obj}{ext}_linesLog_{cycle}.txt'
        results_file = objFolder/f'{obj}{ext}_measurements.txt'
        diagram_file = objFolder/f'{obj}{ext}_diagChart_{cycle}.png'
        print(f'\n- Treating: {obj}{ext}.fits')

        # Load the data
        linesDF = sr.lineslogFile_to_DF(lineLog_file)
        results_dict = sr.loadConfData(results_file, group_variables=False)
        obj_model_conf = obsData[f'{obj}_chemical_model']
        obj_lines_conf = obsData[f'{obj}_line_fitting']

        # Extinction parameters
        # cHbeta = results_dict[f'Extinction_{cycle}']['cHbeta_BR_Hbeta_Hgamma_Hdelta']
        cHbeta = results_dict['Initial_values']['cHbeta_BR_Hbeta_Hgamma_Hdelta']

        # Declare lines to be used in analysis
        idcs_obs = ~linesDF.index.str.contains('_b')
        lineLabels = linesDF.loc[idcs_obs].index
        f_lambda = linesDF.loc[idcs_obs, 'f_lambda'].values
        obsFlux_Err = linesDF.loc[idcs_obs, 'obsFlux':'obsFluxErr'].values
        obsInt_Err = linesDF.loc[idcs_obs, 'obsInt':'obsIntErr'].values

        # Compute matrix fluxes
        cm = Standard_DirectMetchod(n_steps=5000)
        flux_dict = cm.declare_line_fluxes(lineLabels, obsFlux_Err[:, 0], obsFlux_Err[:, 1])

        # Establish extinction correction
        int_dict = cm.red_corr(flux_dict, cHbeta=cHbeta[0], cHbeta_err=cHbeta[1], f_lambda=f_lambda)

        # Plot diagnostics diagram for galaxy
        cm.plot_diag_chart(int_dict, plot_address=diagram_file)

        # Confirm temperature and density diagnostics
        cm.electron_diagnostics(int_dict, neSII_limit_check=obj_model_conf['nSII_lowerDist_check'],
                                Thigh_diag=obj_model_conf['Te_high_diag'])

        # Compute ionic abundances
        cm.ionic_abundances(int_dict, obj_lines_conf, obj_model_conf)

        # Compute elemental abundances


        # Save results abundances
        cm.save_measurments(results_file, cycle_ref)
