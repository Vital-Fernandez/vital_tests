import numpy as np
import pandas as pd
import src.specsiser as sr
import lime

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import chemical_lines_indexing
from astropy.io import fits


def plotting_the_traces(db_address, fits_ext, output_folder, model_ext=''):

    # Load the results
    fit_pickle = sr.load_fit_results(db_address, ext_name=fits_ext, output_format='fits')
    inLines = fit_pickle[f'{fits_ext}_inputs'][0]['line_list']
    inParameters = fit_pickle[f'{fits_ext}_outputs'][0]['parameters_list']
    inFlux = fit_pickle[f'{fits_ext}_inputs'][0]['line_fluxes']
    inErr = fit_pickle[f'{fits_ext}_inputs'][0]['line_err']
    traces_dict = fit_pickle[f'{fits_ext}_traces'][0]

    # Print the results
    print('-- Model parameters table')
    figure_file = f'{output_folder}/{fits_ext}_{model_ext}_fitted_fluxes'
    sr.table_fluxes(figure_file, inLines, inFlux, inErr, traces_dict, merge_dict)

    # Print the results
    print('-- Fitted fluxes table')
    figure_file = f'{output_folder}/{fits_ext}_{model_ext}_MeanOutputs'
    sr.table_params(figure_file, inParameters, traces_dict)

    print('-- Model parameters posterior diagram')
    figure_file = f'{output_folder}/{fits_ext}_{model_ext}_traces_plot.png'
    sr.plot_traces(figure_file, inParameters, traces_dict)

    print('-- Line flux posteriors')
    figure_file = f'{output_folder}/{fits_ext}_{model_ext}_fluxes_grid.png'
    sr.plot_flux_grid(figure_file, inLines, inFlux, inErr, traces_dict)

    print('-- Model parameters corner diagram')
    figure_file = f'{output_folder}/{fits_ext}_{model_ext}_corner.png'
    sr.plot_corner(figure_file, inParameters, traces_dict)

    return


# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

# Direct methods params
merge_dict = {'O2_7319A_b': 'O2_7319A-O2_7330A'}
R_v = obsData['Extinction']['R_v']
red_law = obsData['Extinction']['red_law']

# Grid sampling params
ref_simulations = ['localErr', 'HIICHImistry', 'maxErr',  'noOII']
tech_label = 'GridSampling'

# Load the photoionization grid
file_address_epm = f'{dataFolder}/formated_log_C17_Popstar_1Myr.dat'
grid_epm = pd.read_csv(file_address_epm, delim_whitespace=True, header=0)

gw = sr.GridWrapper()
model_variables = ['logOH', 'logU', 'logNO']
grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_epm, axes_columns=model_variables)
grid_interp = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')

# Pixel to analyse
target_region = [2]
target_voxel = (171, 173)

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = f'{objFolder}/chemistry'
    obsLog_addresss = objFolder / f'{obj}_linesLog.fits'
    absLog_addresss = objFolder / f'{obj}_linesLog_abs.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Loop throught the line regions
    for idx_region in target_region:

        # Voxel mask
        region_label = f'mask_{idx_region}'
        region_mask = fits.getdata(maskFits_address, region_label, ver=1)
        region_mask = region_mask.astype(bool)
        idcs_voxels = np.argwhere(region_mask)

        # Simulation time estimation
        n_voxels = idcs_voxels.shape[0]
        print(f'\nTreating {region_label} consisting of {n_voxels}. The estimated time is {(n_voxels*1.85)/60:0.1f} hrs')

        # Region chemical configuration
        chem_conf_file = dataFolder/f'{obj}_chemical_model_region_{idx_region}.txt'
        chem_conf = lime.load_cfg(chem_conf_file)

        # Compute emissivity grids from the candidate lines
        input_lines = chem_conf['inference_model_configuration']['input_lines_list']
        emis_grid_interp = sr.emissivity_grid_calc(lines_array=input_lines, comp_dict=merge_dict)

        # Confirm in a voxel
        for idx_voxel, idx_pair in enumerate(idcs_voxels):
            idx_j, idx_i = idx_pair
            if idx_j == target_voxel[0] and idx_i == target_voxel[1]:
                print(f'\nTreating voxel {idx_j}-{idx_i}: ({idx_voxel}/{n_voxels})')

                # -------------------------------------  Direct method -------------------------------------------------
                # Output file
                output_direct_method_Db = objFolder / f'tests' / f'{obj}_direct_method_voxel.fits'
                chem_ref = f'{idx_j}-{idx_i}_direct_method'

                # Load voxel data:
                ext_ref = f'{idx_j}-{idx_i}_linelog'
                obs_log = lime.load_lines_log(obsLog_addresss, ext_ref)
                abs_log = lime.load_lines_log(absLog_addresss, ext_ref)

                # Establish and normalize the lines we want
                linesDF = chemical_lines_indexing(input_lines, obs_log, abs_log, obsData)
                lineLabels = linesDF.index.values
                lineWaves = linesDF.wavelength.values
                lineIons = linesDF.ion.values
                lineFluxes = linesDF.obsFlux.values
                lineErr = linesDF.obsFluxErr.values
                lineflambda = sr.flambda_calc(lineWaves, R_V=R_v, red_curve=red_law)

                # Declare object
                obj1_model = sr.SpectraSynthesizer(emis_grid_interp)

                # Declare simulation physical properties
                obj1_model.define_region(lineLabels, lineFluxes, lineErr, lineflambda)

                # Sampling configuration
                obj1_model.simulation_configuration(prior_conf_dict=chem_conf['priors_configuration'],
                                                    highTempIons=chem_conf['simulation_properties']['high_temp_ions_list'])
                # Theoretical model
                obj1_model.inference_model()
                obj1_model.run_sampler(500, 2000, nchains=10, njobs=10)
                obj1_model.save_fit(output_direct_method_Db, ext_name=chem_ref, output_format='fits')

                # Plot the trace
                plotting_the_traces(output_direct_method_Db, chem_ref, objFolder/f'tests', model_ext='direct')

                # # -------------------------------------  Grid sampling -------------------------------------------------
                for j, conf in enumerate(ref_simulations):
                    output_modelFitting_Db = objFolder/'tests'/f'{obj}_{tech_label}_{conf}_voxel.fits'
                    chem_ref = f'{idx_j}-{idx_i}_model_fitting'

                    # Load region fluxes
                    int_DF = lime.load_lines_log(objFolder / f'region_{idx_region}_gridSampling_intensities.txt')
                    err_DF = lime.load_lines_log(objFolder / f'region_{idx_region}_gridSampling_errors.txt')

                    ext_lines = f'{idx_j}-{idx_i}_linelog'
                    int_series = int_DF.loc[ext_lines]
                    err_series = err_DF.loc[ext_lines]
                    region_lines = obsData[f'{tech_label}_{conf}_conf'][f'MASK_{idx_region}_line_list']

                    # Select requested lines, non nan
                    idcs_obs = ~pd.isnull(int_series) & (int_series.index != 'mask') & (int_series.index.isin(region_lines))
                    lineLabels = int_series[idcs_obs].keys().values
                    lineInts = int_series[idcs_obs].values
                    LineErrs = err_series[idcs_obs].values

                    # Error definition according to the model
                    minErr_model = 0.02 if conf != 'maxErr' else np.max(LineErrs / lineInts)

                    # Define model sampler
                    obj1_model = sr.SpectraSynthesizer(grid_sampling=True, grid_interp=grid_interp)
                    obj1_model.define_region(lineLabels, lineInts, LineErrs, minErr=minErr_model)
                    obj1_model.simulation_configuration(prior_conf_dict=obsData['GridSampling_priors'])
                    obj1_model.photoionization_sampling(model_variables)
                    obj1_model.run_sampler(500, 2000, nchains=10, njobs=10, init='advi')
                    obj1_model.save_fit(output_modelFitting_Db, chem_ref, output_format='fits')

                    # Plot the trace
                    plotting_the_traces(output_modelFitting_Db, chem_ref, objFolder/f'tests', model_ext=conf)