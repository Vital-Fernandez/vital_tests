import numpy as np
import pandas as pd
import src.specsiser as sr
import lime

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import chemical_lines_indexing
from astropy.io import fits

import numpy as np
import pyneb as pn

S2 = pn.Atom('S', 2)


# def density_calc(S2, S2_err, S2_atom, temp=15000):
#
#     ratio_array = np.random.normal(S2, S2_err, 1000)
#     ne_array = S2_atom.getTemDen(ratio_array, tem=temp, to_eval='L(6717)/L(6731)')
#     print(f'n_e = {np.nanmean(ne_array):.2f}+/-{np.nanstd(ne_array):.2f}')
#
#     return
#
# density_calc(1.40, 0.1, S2)
#
# density_calc(1.40, 0.1, S2)
# density_calc(1.31, 0.09, S2)
# density_calc(1.34, 0.1, S2)
#0.134  0.003 0.200  0.006
#0.696  0.026 0.538  0.052
rc = pn.RedCorr(R_V=3.1, E_BV=0.15, law='G03 LMC')

S2_array = np.random.normal(0.200, 0.006, 1000) * 1e-6
S3_array = np.random.normal(0.538, 0.052, 1000) * 1e-6
O_log = 12 + np.log10(S2_array+S3_array)
print(f'S/H = {O_log.mean():0.2f} +/- {O_log.std():0.2f}')

def Ebv_from_cHbeta(self, cHbeta, reddening_curve, R_v):
    if cHbeta == None:
        exit('Warning: no cHbeta or E(B-V) provided to reddening curve, code aborted')
    else:
        if cHbeta != None:
            E_BV = cHbeta * 2.5 / self.reddening_Xx(array([self.Hbeta_wavelength]), reddening_curve, R_v)[0]
            return E_BV


def flambda_from_Xx(self, Xx, reddening_curve, R_v):
    X_Hbeta = self.reddening_Xx(array([self.Hbeta_wavelength]), reddening_curve, R_v)[0]

    f_lines = Xx / X_Hbeta - 1

    return f_lines


def add_abundances(O2, O2_err, O3, O3_err, size=1000):

    O2_array = np.random.normal(loc=O2, scale=O2_err, size=size)
    O3_array = np.random.normal(loc=O3, scale=O3_err, size=size)

    O_array_nat = np.power(10, O2_array-12) + np.power(10, O3_array-12)
    O_array = 12 + np.log10(O_array_nat)

    print(f'12 + log(O/H) = {O_array.mean():.2f} +/- {O_array.std():.2f}')

    return


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
ref_simulations = ['localErr', 'HIICHImistry', 'maxErr',  'minOneErr']
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
target_voxel = (147, 162)

order_array = np.array([0.1, 1, 10, 100, 1000]) / 100

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

                # # -------------------------------------  Direct method -------------------------------------------------
                # # Output file
                # output_direct_method_Db = objFolder / f'tests' / f'{obj}_direct_method_voxel.fits'
                # chem_ref = f'{idx_j}-{idx_i}_direct_method'
                #
                # # Load voxel data:
                # ext_ref = f'{idx_j}-{idx_i}_linelog'
                # obs_log = lime.load_lines_log(obsLog_addresss, ext_ref)
                # abs_log = lime.load_lines_log(absLog_addresss, ext_ref)
                #
                # # Establish and normalize the lines we want
                # linesDF = chemical_lines_indexing(input_lines, obs_log, abs_log, obsData)
                # lineLabels = linesDF.index.values
                # lineWaves = linesDF.wavelength.values
                # lineIons = linesDF.ion.values
                # lineFluxes = linesDF.obsFlux.values
                # lineErr = linesDF.obsFluxErr.values
                # lineflambda = sr.flambda_calc(lineWaves, R_V=R_v, red_curve=red_law)
                #
                # # Declare object
                # obj1_model = sr.SpectraSynthesizer(emis_grid_interp)
                #
                # # Declare simulation physical properties
                # obj1_model.define_region(lineLabels, lineFluxes, lineErr, lineflambda)
                #
                # # Sampling configuration
                # obj1_model.simulation_configuration(prior_conf_dict=chem_conf['priors_configuration'],
                #                                     highTempIons=chem_conf['simulation_properties']['high_temp_ions_list'])
                # # Theoretical model
                # obj1_model.inference_model()
                # obj1_model.run_sampler(500, 2000, nchains=10, njobs=10)
                # obj1_model.save_fit(output_direct_method_Db, ext_name=chem_ref, output_format='fits')
                #
                # # Plot the trace
                # plotting_the_traces(output_direct_method_Db, chem_ref, objFolder/f'tests', model_ext='direct')

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
                    # minErr_model = 0.02 if conf not in ['maxErr', 'joinOIIMax'] else np.max(LineErrs / lineInts)

                    if conf == 'maxErr':
                        minErr_in = np.max(LineErrs / lineInts)
                    elif conf == 'minOneErr':
                        max_err = np.max(LineErrs / lineInts)
                        order_array = np.array([0.002, 0.02, 0.2, 2, 10, 100, 1000, 10000]) / 100
                        idx_order = np.argmin(np.abs(order_array - max_err))
                        min_err = order_array[idx_order - 1]
                        minErr_model = LineErrs / lineInts

                        idcs_smaller_err = minErr_model < min_err
                        minErr_model[idcs_smaller_err] = min_err
                        LineErrs = minErr_model * lineInts
                        minErr_in = None
                    else:
                        minErr_in = 0.02

                    # Define model sampler
                    obj1_model = sr.SpectraSynthesizer(grid_sampling=True, grid_interp=grid_interp)
                    obj1_model.define_region(lineLabels, lineInts, LineErrs, minErr=minErr_in)
                    obj1_model.simulation_configuration(prior_conf_dict=obsData['GridSampling_priors'])
                    obj1_model.photoionization_sampling(model_variables)
                    obj1_model.run_sampler(500, 2000, nchains=10, njobs=10, init='advi')
                    obj1_model.save_fit(output_modelFitting_Db, chem_ref, output_format='fits')

                    # Plot the trace
                    plotting_the_traces(output_modelFitting_Db, chem_ref, objFolder/f'tests', model_ext=conf)