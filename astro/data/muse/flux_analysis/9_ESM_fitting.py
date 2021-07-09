import numpy as np
import pandas as pd
from pathlib import Path
import src.specsiser as sr
from src.specsiser.inference_model import fits_db
from astro.data.muse.common_methods import grid_columns
from astropy.io import fits

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

# Load the data
grid_file = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/HII-CHI-mistry_1Myr_grid.csv')
grid_DF = pd.read_csv(grid_file, skiprows=1, names=grid_columns.values())
grid_DF.logNO = np.round(grid_DF.logNO.values, decimals=3)
df_columns = np.array(['wavelength', 'ion',	'intg_flux', 'intg_err'])

# Remove carbon dimension
idcs_rows = grid_DF.carbon == 'O'
idcs_columns = ~grid_DF.columns.isin(['carbon'])
grid_3D_DF = grid_DF.loc[idcs_rows, idcs_columns]

model_variables = ['logOH', 'logU', 'logNO']
gw = sr.ModelGridWrapper()
grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_3D_DF, axes_columns=model_variables)
grid_interpolators = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')

idx = 1
print('logOH', axes_cords_a['logOH'][idx])
print('logU', axes_cords_a['logU'][idx])
print('logNO', axes_cords_a['logNO'][idx])
print('Dictionary grid', grid_dict['O2_3729A'][idx, idx, idx])

exclude_lines = np.array(['H1_4861A', 'N1_5198A', 'N1_5200A', 'C2_4267A', 'O2_4651A', 'C2_4659A'])
min_wavelength = 4700
input_lines = np.array(list(grid_dict.keys()))

ion_array, wave_array, latex_array = sr.label_decomposition(input_lines)
idcs_lines = ~np.isin(input_lines, exclude_lines) & (min_wavelength < wave_array)
input_lines = input_lines[idcs_lines]

# cord_true = {'logOH': 7.2, 'logU': -3.75, 'logNO': -1.875}
cord_true = {'logOH': 8.255, 'logU': -2.155, 'logNO': -1.950}

logOH_range = np.array([cord_true['logOH']])
logU_range = np.array([cord_true['logU']])
logNO_range = np.array([cord_true['logNO']])

# Default configuration
conf_file = Path(dataFolder/'CGCG007/chemical_fitting_conf.txt')
conf_fit = sr.loadConfData(conf_file)

# Loop throught the grid of synthetic conditions
for i, logOH in enumerate(logOH_range):
    for j, logU in enumerate(logU_range):
        for k, logNO in enumerate(logNO_range):

            # True value coordinate for interpolation
            coord_true = [[logOH, logU, logNO]]
            header_params = {'logOH': logOH, 'logU': logU, 'logNO': logNO}

            # Output files
            objFolder = resultsFolder/f'CGCG007/'
            cord_label = f'{logOH*1000:.0f}_{logU*-1000:.0f}_{logNO*-1000:.0f}'
            outputCfg = objFolder/f'{cord_label}.txt'
            outputDb = objFolder/f'{cord_label}'
            outputFits = objFolder/f'{cord_label}.fits'

            # Fill the dataframe with integrated flux
            linesDF = pd.DataFrame(index=input_lines, columns=df_columns)
            for line in input_lines:
                ion, wavelength, latexLabel = sr.label_decomposition(line, scalar_output=True)
                flux = np.power(10, grid_interpolators[line](coord_true).eval()[0][0])
                linesDF.loc[line, :] = wavelength, ion, flux, flux * 0.05

            # Declare sampler
            obj1_model = sr.SpectraSynthesizer()

            # Declare region physical model
            obj1_model.define_region(linesDF)

            # Declare region physical model
            obj1_model.simulation_configuration(model_parameters=conf_fit['inference_model_configuration']['parameter_list'],
                                                prior_conf_dict=conf_fit['priors_configuration'],
                                                grid_interpolator=grid_interpolators)

            obj1_model.photoionization_sampling(conf_fit['inference_model_configuration']['parameter_list'])

            obj1_model.run_sampler(3000, 2000, nchains=3, njobs=3)

            obj1_model.save_fit(outputCfg, cord_label, output_format='cfg')
            obj1_model.save_fit(outputDb, cord_label, output_format='pickle')
            obj1_model.save_fit(outputFits, cord_label, output_format='fits', user_header=header_params)

            # Plot the results

            # fit_pickle = sr.load_MC_fitting(outputDb, output_format='pickle')
            # inLines, inParameters = fit_pickle['inputs']['line_list'], fit_pickle['inputs']['parameter_list']
            # inFlux, inErr = fit_pickle['inputs']['line_fluxes'].astype(float), fit_pickle['inputs']['line_err'].astype(float)
            # traces_dict = fit_pickle['outputs']

            # Load the results
            fit_fits = sr.load_MC_fitting(outputFits, output_format='fits')
            inLines, inParameters = fit_fits['inputs'][0]['line_list'], fit_fits['outputs'][0]['parameters_list']
            inFlux, inErr = fit_fits['inputs'][0]['line_fluxes'], fit_fits['inputs'][0]['line_err']
            traces_dict = fit_fits['traces'][0]

            # print('-- Model parameters table')
            figure_file = objFolder / f'{cord_label}_pI_fitting_MeanOutputs'
            obj1_model.table_mean_outputs(figure_file, inParameters, traces_dict, true_values=header_params)

            print('-- Model parameters posterior diagram')
            figure_file = objFolder / f'{cord_label}_pI_fitting_ParamsPosteriors.png'
            obj1_model.tracesPosteriorPlot(figure_file, inParameters, traces_dict)

            print('-- Model parameters corner diagram')
            figure_file = objFolder / f'{cord_label}_pI_fitting_cornerPlot.png'
            obj1_model.corner_plot(figure_file, inParameters, traces_dict, true_values=header_params)

            print('-- Model emission flux posteriors')
            figure_file = objFolder/f'{cord_label}_pI_EmFluxPosteriors.png'
            obj1_model.fluxes_distribution(figure_file, inLines, inFlux, inErr, traces_dict)

            print('-- Model emission flux table')
            figure_file = objFolder/f'{cord_label}_pI_EmFluxPosteriors'
            obj1_model.table_line_fluxes(figure_file, inLines, inFlux, inErr, traces_dict)
