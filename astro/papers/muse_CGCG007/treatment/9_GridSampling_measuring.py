import numpy as np
import pandas as pd
import src.specsiser as sr
import lime

from pathlib import Path
from astropy.io import fits

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

# Measurement files and reference files
ref_simulations = ['localErr', 'HIICHImistry', 'maxErr',  'noOII']
tech_label = 'GridSampling'

# Load the photoionization grid
model_variables = ['logOH', 'logU', 'logNO']

file_address = f'{dataFolder}/HII-CHI-mistry_1Myr_grid_O.txt'
grid_3D = pd.read_csv(file_address, delim_whitespace=True, header=0)
# gw = sr.GridWrapper()
# grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_3D, axes_columns=model_variables)
# grid_interp = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')

file_address_epm = f'{dataFolder}/formated_log_C17_Popstar_1Myr.dat'
grid_epm = pd.read_csv(file_address_epm, delim_whitespace=True, header=0)
gw = sr.GridWrapper()
grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_epm, axes_columns=model_variables)
grid_interp = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')

# Loop throught the objects and masks
for i, obj in enumerate(objList):

    # Input data
    objFolder = resultsFolder/obj
    db_addresss = objFolder / f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    chemFolder = objFolder / 'chemistry'

    for j, conf in enumerate(ref_simulations):

        outputFits = objFolder / f'NewGrid_{obj}_{tech_label}_{conf}.fits'

        # Loop throught the line regions
        for idx_region in [0, 1, 2]:

            # Input lines
            region_lines = obsData[f'{tech_label}_{conf}_conf'][f'MASK_{idx_region}_line_list']

            # Load voxel list
            region_label = f'mask_{idx_region}'
            region_mask = fits.getdata(maskFits_address, region_label, ver=1)
            region_mask = region_mask.astype(bool)
            idcs_voxels = np.argwhere(region_mask)
            n_voxels = idcs_voxels.shape[0]

            # Load region fluxes
            int_DF = lime.load_lines_log(objFolder/f'region_{idx_region}_gridSampling_intensities.txt')
            err_DF = lime.load_lines_log(objFolder/f'region_{idx_region}_gridSampling_errors.txt')

            # Loop through the region voxels
            for idx_voxel, idx_pair in enumerate(idcs_voxels):

                idx_j, idx_i = idx_pair
                ext_lines = f'{idx_j}-{idx_i}_linelog'
                ext_chem = f'{idx_j}-{idx_i}_gridSampler'
                print(f'\nTreating voxel {idx_j}-{idx_i}: ({tech_label}_{conf}: ({idx_voxel}/{n_voxels})\n')

                # Load voxel fluxes:
                int_series = int_DF.loc[ext_lines]
                err_series = err_DF.loc[ext_lines]

                # Select requested lines, non nan
                idcs_obs = ~pd.isnull(int_series) & (int_series.index != 'mask') & (int_series.index.isin(region_lines))
                lineLabels = int_series[idcs_obs].keys().values
                lineInts = int_series[idcs_obs].values
                LineErrs = err_series[idcs_obs].values

                # Error definition according to the model
                minErr_model = 0.02 if conf != 'maxErr' else np.max(LineErrs/lineInts)

                # Define model sampler
                obj1_model = sr.SpectraSynthesizer(grid_sampling=True, grid_interp=grid_interp)
                obj1_model.define_region(lineLabels, lineInts, LineErrs, minErr=minErr_model)
                obj1_model.simulation_configuration(prior_conf_dict=obsData['GridSampling_priors'])
                obj1_model.photoionization_sampling(model_variables)
                obj1_model.run_sampler(500, 2000, nchains=10, njobs=10, init='advi')
                obj1_model.save_fit(outputFits, ext_chem, output_format='fits')

            # # Load the results
            # fit_results = sr.load_fit_results(outputFits, ext_name=ext_chem, output_format='fits')
            # inLines = fit_results[f'{ext_chem}_inputs'][0]['line_list']
            # inParameters = fit_results[f'{ext_chem}_outputs'][0]['parameters_list']
            # inFlux = fit_results[f'{ext_chem}_inputs'][0]['line_fluxes']
            # inErr = fit_results[f'{ext_chem}_inputs'][0]['line_err']
            # traces_dict = fit_results[f'{ext_chem}_traces'][0]

            # # Print the results
            # print('-- Model parameters table')
            # figure_file = f'{chemFolder}/{ext_chem}_fitted_fluxes'
            # sr.table_fluxes(figure_file, inLines, inFlux, inErr, traces_dict)
            #
            # # Print the results
            # print('-- Fitted fluxes table')
            # figure_file = f'{chemFolder}/{ext_chem}_MeanOutputs'
            # sr.table_params(figure_file, inParameters, traces_dict)
            #
            # print('-- Model parameters posterior diagram')
            # figure_file = f'{chemFolder}/{ext_chem}_traces_plot.png'
            # sr.plot_traces(figure_file, inParameters, traces_dict)
            #
            # print('-- Line flux posteriors')
            # figure_file = f'{chemFolder}/{ext_chem}_fluxes_grid.png'
            # sr.plot_flux_grid(figure_file, inLines, inFlux, inErr, traces_dict)

            # print('-- Model parameters posterior diagram')
            # figure_file = f'{chemFolder}/{ext_chem}_trace_plot.png'
            # sr.plot_traces(figure_file, inParameters, traces_dict)

            # # Load the results
            # ext_chem = f'{idx_j}-{idx_i}_chemistry'
            # fit_results = sr.load_fit_results(outputDb, ext_name=ext_chem, output_format='fits')
            # cHBeta = fit_results[f'{ext_chem}_outputs'][1]['cHbeta']
            # print(f'- cHBeta: {cHBeta}')

            # inLines = fit_results[f'{ext_chem}_inputs'][0]['line_list']
            # inParameters = fit_results[f'{ext_chem}_outputs'][0]['parameters_list']
            # inFlux = fit_results[f'{ext_chem}_inputs'][0]['line_fluxes']
            # inErr = fit_results[f'{ext_chem}_inputs'][0]['line_err']
            # traces_dict = fit_results[f'{ext_chem}_traces'][0]
