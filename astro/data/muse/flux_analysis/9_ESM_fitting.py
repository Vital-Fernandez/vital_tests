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

Rv = obsData['Extinction']['R_v']
red_law = obsData['Extinction']['red_law']

# # Load the data
# grid_file = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/HII-CHI-mistry_1Myr_grid.csv')
# grid_DF = pd.read_csv(grid_file, skiprows=1, names=grid_columns.values())
#
# grid_DF.logNO = np.round(grid_DF.logNO.values, decimals=3)
# df_columns = np.array(['wavelength', 'ion',	'intg_flux', 'intg_err'])
#
# for column in grid_DF.columns:
#     if column in grid_columns:
#         print(f'{column} -> {grid_columns[column]}')
# # Remove carbon dimension
# idcs_rows = grid_DF.carbon == 'O'
# idcs_columns = ~grid_DF.columns.isin(['carbon'])
# grid_3D_DF = grid_DF.loc[idcs_rows, idcs_columns]
#
# file_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/HII-CHI-mistry_1Myr_grid_O.txt')
# with open(file_address, 'wb') as output_file:
#     string_DF = grid_3D_DF.to_string(index=False)
#     output_file.write(string_DF.encode('UTF-8'))
#
# idx = 1
# print('logOH', axes_cords_a['logOH'][idx])
# print('logU', axes_cords_a['logU'][idx])
# print('logNO', axes_cords_a['logNO'][idx])
# print('Dictionary grid', grid_dict['O2_3729A'][idx, idx, idx])

file_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/HII-CHI-mistry_1Myr_grid_O.txt')
grid_3D_DF = pd.read_csv(file_address, delim_whitespace=True, header=0)

# # Adding line ratios:
# R_S2S3 = (np.power(10, grid_3D_DF['S2_6716A']) + np.power(10, grid_3D_DF['S2_6731A']))\
#          /(np.power(10, grid_3D_DF['S3_9069A']) + np.power(10, grid_3D_DF['S3_9531A']))
# grid_3D_DF['R_S2S3'] = np.log10(R_S2S3)
#
# R_O3N2 = np.power(10, grid_3D_DF['O3_5007A']) / np.power(10, grid_3D_DF['N2_6584A'])
# grid_3D_DF['R_O3N2'] = np.log10(R_O3N2)
# file_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/HII-CHI-mistry_1Myr_grid_O.txt')
#
# with open(file_address, 'wb') as output_file:
#     string_DF = grid_3D_DF.to_string(index=False)
#     output_file.write(string_DF.encode('UTF-8'))

model_variables = ['logOH', 'logU', 'logNO']
gw = sr.ModelGridWrapper()
grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_3D_DF, axes_columns=model_variables)
grid_interpolators = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')

model_lines = np.array(list(grid_dict.keys()))

lines_MUSE = np.array(['O3_4959A', 'O3_5007A',
                       'He1_5876A',
                       'N2_6548A', 'H1_6563A', 'N2_6584A',
                       'S2_6716A', 'S2_6731A',
                       'S3_9069A', 'S3_9531A'])

ratio_format = {'R_S2S3': ['-', '-', '$R_{S2S3}$'],
                'R_O3N2': ['-', '-', '$R_{O3N2}$']}


idcs_lines = np.isin(lines_MUSE, model_lines)
input_lines = lines_MUSE[idcs_lines]

cord_true = {'logOH': 7.2, 'logU': -3.75, 'logNO': -1.875, 'cHbeta': 0.255}
# cord_true = {'logOH': 8.255, 'logU': -2.155, 'logNO': -1.950, 'cHbeta': 0.255}

logOH_range = np.array([cord_true['logOH']])
logU_range = np.array([cord_true['logU']])
logNO_range = np.array([cord_true['logNO']])

# Default configuration
conf_file = Path(dataFolder/'CGCG007/chemical_fitting_conf.txt')
conf_fit = sr.loadConfData(conf_file)

# Extinction model
objRed = sr.ExtinctionModel(Rv=Rv, red_curve=red_law)

# Loop throught the grid of synthetic conditions
for i, logOH in enumerate(logOH_range):
    for j, logU in enumerate(logU_range):
        for k, logNO in enumerate(logNO_range):

            # True value coordinate for interpolation
            coord_true = [[logOH, logU, logNO]]
            header_params = {'logOH': logOH, 'logU': logU, 'logNO': logNO, 'cHbeta': cord_true['cHbeta']}

            # Output files
            objFolder = resultsFolder/f'CGCG007/'
            cord_label = f'{logOH*1000:.0f}_{logU*-1000:.0f}_{logNO*-1000:.0f}'
            outputCfg = objFolder/f'{cord_label}.txt'
            outputDb = objFolder/f'{cord_label}'
            outputFits = objFolder/f'{cord_label}.fits'

            # Fill the dataframe with integrated flux
            linesDF = pd.DataFrame(index=input_lines, columns=['wavelength', 'ion',	'intg_flux', 'intg_err'])
            for line in input_lines:
                if not line.startswith('R_'):
                    ion, wavelength, latexLabel = sr.label_decomposition(line, scalar_output=True)
                    flux = np.power(10, grid_interpolators[line](coord_true).eval()[0][0])
                    f_lambda = objRed.gasExtincParams(wave=wavelength)
                    flux_red = flux * np.power(10, -cord_true['cHbeta'] * f_lambda)
                    linesDF.loc[line, :] = wavelength, ion, flux_red, flux_red * 0.05
                else:
                    ion, wavelength, latexLabel = 'None', 0.0, line
                    flux = np.power(10, grid_interpolators[line](coord_true).eval()[0][0])
                    linesDF.loc[line, :] = wavelength, ion, flux, flux * 0.05

            # Declare sampler
            obj1_model = sr.SpectraSynthesizer()

            # Declare region physical model
            lineLabels = linesDF.index.values
            lineFluxes = linesDF.intg_flux.values
            lineErr = linesDF.intg_err.values
            obj1_model.define_region(lineLabels, lineFluxes, lineErr, extinction_model=objRed)

            # Declare region physical model
            obj1_model.simulation_configuration(model_parameters=conf_fit['inference_model_configuration']['parameter_list'],
                                                prior_conf_dict=conf_fit['priors_configuration'],
                                                grid_interpolator=grid_interpolators)
            obj1_model.photoionization_sampling(conf_fit['inference_model_configuration']['parameter_list'])

            obj1_model.run_sampler(1000, 3000, nchains=8, njobs=4, init='advi')

            # obj1_model.save_fit(outputCfg, cord_label, output_format='cfg')
            obj1_model.save_fit(outputDb, cord_label, output_format='pickle')
            # obj1_model.save_fit(outputFits, cord_label, output_format='fits', user_header=header_params)

            # Plot the results
            fit_pickle = sr.load_MC_fitting(outputDb, output_format='pickle')
            inLines, inParameters = fit_pickle['inputs']['line_list'], fit_pickle['inputs']['parameter_list']
            inFlux, inErr = fit_pickle['inputs']['line_fluxes'].astype(float), fit_pickle['inputs']['line_err'].astype(float)
            traces_dict = fit_pickle['outputs']

            # Load the results
            # fit_fits = sr.load_MC_fitting(outputFits, output_format='fits')
            # in_header, out_header, trace_header = fit_fits['inputs'][1], fit_fits['outputs'][1], fit_fits['traces'][1]
            # inLines, inParameters = fit_fits['inputs'][0]['line_list'], fit_fits['outputs'][0]['parameters_list']
            # inFlux, inErr = fit_fits['inputs'][0]['line_fluxes'], fit_fits['inputs'][0]['line_err']
            # traces_dict = fit_fits['traces'][0]

            # print('-- Model parameters table')
            figure_file = objFolder / f'{cord_label}_pI_fitting_MeanOutputs'
            obj1_model.table_mean_outputs(figure_file, inParameters, traces_dict, true_values=header_params)

            print('-- Model parameters posterior diagram')
            figure_file = objFolder / f'{cord_label}_pI_fitting_ParamsPosteriors.png'
            obj1_model.tracesPosteriorPlot(figure_file, inParameters, traces_dict, true_values=header_params)

            print('-- Model emission flux posteriors')
            figure_file = objFolder/f'{cord_label}_pI_EmFluxPosteriors.png'
            obj1_model.fluxes_distribution(figure_file, inLines, inFlux, inErr, traces_dict, user_labels=ratio_format)

            print('-- Model emission flux table')
            figure_file = objFolder/f'{cord_label}_pI_EmFluxPosteriors'
            obj1_model.table_line_fluxes(figure_file, inLines, inFlux, inErr, traces_dict, user_labels=ratio_format)

            print('-- Model parameters corner diagram')
            figure_file = objFolder / f'{cord_label}_pI_fitting_cornerPlot.png'
            obj1_model.corner_plot(figure_file, inParameters, traces_dict, true_values=header_params)