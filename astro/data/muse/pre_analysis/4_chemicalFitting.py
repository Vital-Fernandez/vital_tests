import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, red_corr_HalphaHbeta_ratio
from timeit import default_timer as timer

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['sample_data']['object_list']
fileList = obsData['sample_data']['file_list']
model_conf_file = obsData['sample_data']['file_list']
dataFolder = Path(obsData['sample_data']['data_folder'])
resultsFolder = Path(obsData['sample_data']['results_folder'])
z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']
dict_errs = {}

start = timer()
for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        objFolder = resultsFolder/obj
        voxelFolder = resultsFolder/obj/'voxel_data'
        db_addresss = objFolder/f'{obj}_database.txt'
        cube_address = dataFolder/fileList[i]
        mask_address = objFolder/f'{obj}_mask.txt'
        chem_conf_file = resultsFolder/obj/'chemical_model_config.txt'
        chem_conf = sr.loadConfData(chem_conf_file, group_variables=False)

        # Load data
        obj_db = pd.read_csv(db_addresss, delim_whitespace=True, header=0, index_col=0)

        # Add columns in the object database for the chemical fitting results
        for param in chem_conf['inference_model_configuration']['parameter_list']:
            obj_db[param] = np.nan
            obj_db[param + '_err'] = np.nan

        wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')

        # Get line region data
        lineFlux_dict, levelFlux_dict, levelText_dict = compute_line_flux_image(lineAreas,
                                                                                cube,
                                                                                z_objs[i],
                                                                                percent_array=pertil_array)

        # Declare voxels to analyse
        fluxImage_6312 = lineFlux_dict['S3_6312A']
        fluxLevels_6312 = levelFlux_dict['S3_6312A']
        int_level = fluxLevels_6312[-4]#[-4]
        array_idcs_voxels = np.argwhere(fluxImage_6312 > int_level)

        for idx_voxel, idx_pair in enumerate(array_idcs_voxels):

            idx_j, idx_i = idx_pair

            # Data location
            outputDb = voxelFolder/f'{idx_j}-{idx_i}_fitting.db'
            local_lineslog = voxelFolder/f'{idx_j}-{idx_i}_lineslog_{obj}.txt'
            inputLinesLog = voxelFolder/f'{idx_j}-{idx_i}_inputLinesLog.txt'

            # Asociate the corresponding flux to the appropiate line
            linesDF = sr.lineslogFile_to_DF(local_lineslog)

            if 'H1_4861A' in linesDF.index:

                # Sampler object
                obj1_model = sr.SpectraSynthesizer()

                # Adjust the lines log for the fitting
                flux_Hbeta = linesDF.loc['H1_4861A', 'intg_flux']
                linesDF.insert(loc=1, column='obsFlux', value=np.nan)
                linesDF.insert(loc=2, column='obsFluxErr', value=np.nan)
                linesDF['obsFlux'] = linesDF['intg_flux']/flux_Hbeta
                linesDF['obsFluxErr'] = linesDF['intg_err']/flux_Hbeta
                linesDF.rename(columns={'wavelength': 'obsWave'}, inplace=True)
                for ion in ['H1', 'He1', 'He2']:
                    idcs_ions = linesDF.ion == ion
                    linesDF.loc[idcs_ions, 'ion'] = ion + 'r'
                if ('O2_7319A' in linesDF.index) and ('O2_7330A' in linesDF.index):
                    flux_comb = linesDF.loc['O2_7319A', 'intg_flux'] + linesDF.loc['O2_7330A', 'intg_flux']
                    err_comb = np.sqrt(linesDF.loc['O2_7319A', 'intg_err']**2 + linesDF.loc['O2_7330A', 'intg_err']**2)
                    linesDF.loc['O2_7319A_b'] = None
                    linesDF.loc['O2_7319A_b', ['obsWave', 'obsFlux', 'obsFluxErr']] = 7325, flux_comb/flux_Hbeta, err_comb/flux_Hbeta
                    linesDF.loc['O2_7319A_b', 'ion'] = 'O2'
                sr.save_lineslog(linesDF, inputLinesLog)

                # Load emission lines
                objLinesDF = sr.import_emission_line_data(inputLinesLog,
                                                    include_lines=chem_conf['inference_model_configuration']['input_lines'])



                # Declare simulation physical properties
                objRed = sr.ExtinctionModel(Rv=chem_conf['simulation_properties']['R_v'],
                                            red_curve=chem_conf['simulation_properties']['reddenig_curve'],
                                            data_folder=chem_conf['data_location']['external_data_folder'])

                objIons = sr.IonEmissivity(tempGrid=chem_conf['simulation_properties']['temp_grid'],
                                           denGrid=chem_conf['simulation_properties']['den_grid'])

                # Generate interpolator from the emissivity grids
                ionDict = objIons.get_ions_dict(np.unique(objLinesDF.ion.values))
                objIons.computeEmissivityGrids(objLinesDF, ionDict, linesDb=sr._linesDb,
                                                combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})

                # Declare chemical model
                objChem = sr.DirectMethod(linesDF=objLinesDF,
                                          highTempIons=chem_conf['simulation_properties']['high_temp_ions_list'])

                # Declare region physical model
                obj1_model.define_region(objLinesDF, objIons, objRed, objChem)

                # Declare sampling properties
                obj1_model.simulation_configuration(chem_conf['inference_model_configuration']['parameter_list'],
                                                    prior_conf_dict=chem_conf['priors_configuration'])

                # Declare simulation inference model
                obj1_model.inference_model(include_Thigh_prior=chem_conf['inference_model_configuration']['T_high_check'])

                # Run the simulation
                obj1_model.run_sampler(outputDb, 5000, 2000, njobs=1)

                if outputDb.is_file() and Path(outputDb).with_suffix('.txt').is_file():

                    # Plot the results
                    fit_results = sr.load_MC_fitting(outputDb)

                    # Print the results
                    print('-- Model parameters table')
                    figure_file = voxelFolder / f'{idx_j}-{idx_i}_MeanOutputs'
                    obj1_model.table_mean_outputs(figure_file, fit_results)

                    print('-- Flux values table')
                    figure_file = voxelFolder / f'{idx_j}-{idx_i}_FluxComparison'
                    obj1_model.table_line_fluxes(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})

                    print('-- Model parameters posterior diagram')
                    figure_file = voxelFolder / f'{idx_j}-{idx_i}_ParamsPosteriors.png'
                    obj1_model.tracesPosteriorPlot(figure_file, fit_results)

                    print('-- Line flux posteriors')
                    figure_file = voxelFolder / f'{idx_j}-{idx_i}_lineFluxPosteriors.png'
                    obj1_model.fluxes_distribution(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})
end = timer()
print(f'Working time:{end-start}')