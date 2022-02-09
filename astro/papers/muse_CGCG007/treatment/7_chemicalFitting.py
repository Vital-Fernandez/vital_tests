import numpy as np
import src.specsiser as sr
import time
import lime

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import voxel_security_check
from astropy.io import fits

import os
# Declare data and files location
obsData = sr.loadConfData('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

merge_dict = {'O2_7319A_b': 'O2_7319A-O2_7330A'}

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    db_addresss = objFolder / f'{obj}_database.fits'
    fitsLog_addresss = objFolder / f'{obj}_linesLog.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Output files
    outputDb = objFolder/f'{obj}_chemical.fits'

    # Loop throught the line regions
    start = time.time()
    for idx_region in [0, 1, 2]:
    # for idx_region in [0]:

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
        chem_conf = sr.loadConfData(chem_conf_file, group_variables=False)

        # Load emission lines
        input_lines = chem_conf['inference_model_configuration']['input_lines']
        ion_array, wave_array, latex_array = lime.label_decomposition(input_lines)

        objIons = sr.IonEmissivity(tempGrid=chem_conf['simulation_properties']['temp_grid'],
                                   denGrid=chem_conf['simulation_properties']['den_grid'])

        ionDict = objIons.get_ions_dict(ion_array)

        # Generate interpolator from the emissivity grids
        objIons.computeEmissivityGrids(input_lines, ionDict, combined_dict=merge_dict)

        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            idx_j, idx_i = idx_pair
            print(f'\nTreating voxel {idx_j}-{idx_i}: ({idx_voxel}/{n_voxels})')

            # Data location
            chem_ref = f'{idx_j}-{idx_i}_chemistry'

            # Load voxel data:
            ext_ref = f'{idx_j}-{idx_i}_linelog'
            linesDF = lime.load_lines_log(fitsLog_addresss, ext_ref)

            if voxel_security_check(linesDF):

                linesDF.insert(loc=1, column='obsFlux', value=np.nan)
                linesDF.insert(loc=2, column='obsFluxErr', value=np.nan)

                idcs_gaussian = (linesDF.blended_label != 'None') & (~linesDF.index.str.contains('_b'))
                linesDF.loc[~idcs_gaussian, 'obsFlux'] = linesDF.loc[~idcs_gaussian, 'intg_flux']
                linesDF.loc[~idcs_gaussian, 'obsFluxErr'] = linesDF.loc[~idcs_gaussian, 'intg_err']
                linesDF.loc[idcs_gaussian, 'obsFlux'] = linesDF.loc[idcs_gaussian, 'gauss_flux']
                linesDF.loc[idcs_gaussian, 'obsFluxErr'] = linesDF.loc[idcs_gaussian, 'gauss_err']

                # Normalize by Hbeta the lines log for the fitting
                flux_Hbeta = linesDF.loc['H1_4861A', 'obsFlux']
                linesDF['obsFlux'] = linesDF['obsFlux'] / flux_Hbeta
                linesDF['obsFluxErr'] = linesDF['obsFluxErr'] / flux_Hbeta

                # linesDF.rename(columns={'wavelength': 'obsWave'}, inplace=True)
                if ('O2_7319A' in linesDF.index) and ('O2_7330A' in linesDF.index):
                    flux_comb = linesDF.loc['O2_7319A', 'obsFlux'] + linesDF.loc['O2_7330A', 'obsFlux']
                    err_comb = np.sqrt(linesDF.loc['O2_7319A', 'obsFluxErr']**2 + linesDF.loc['O2_7330A', 'obsFluxErr']**2)
                    linesDF.loc['O2_7319A_b'] = None
                    linesDF.loc['O2_7319A_b', ['wavelength', 'obsFlux', 'obsFluxErr']] = 7325, flux_comb, err_comb
                    linesDF.loc['O2_7319A_b', 'ion'] = 'O2'

                objLinesDF = sr.import_emission_line_data(linesDF=linesDF, include_lines=input_lines)

                normLine = 'H1_4861A'
                idcs_lines = (objLinesDF.index != normLine)
                lineLabels = objLinesDF.loc[idcs_lines].index
                lineIons = objLinesDF.loc[idcs_lines, 'ion'].values
                lineFluxes = objLinesDF.loc[idcs_lines, 'obsFlux'].values
                lineErr = objLinesDF.loc[idcs_lines, 'obsFluxErr'].values

                # Declare simulation physical properties
                objRed = sr.ExtinctionModel(Rv=chem_conf['simulation_properties']['R_v'],
                                            red_curve=chem_conf['simulation_properties']['reddenig_curve'],
                                            data_folder=chem_conf['data_location']['external_data_folder'])

                # Declare chemical model
                objChem = sr.DirectMethod(lineLabels, highTempIons=chem_conf['simulation_properties']['high_temp_ions_list'])

                # Prepare simulation
                obj1_model = sr.SpectraSynthesizer()

                obj1_model.define_region(lineLabels, lineFluxes, lineErr, objIons, objRed, objChem)

                obj1_model.simulation_configuration(chem_conf['inference_model_configuration']['parameter_list'],
                                                    prior_conf_dict=chem_conf['priors_configuration'],
                                                    photo_ionization_grid=False)

                obj1_model.simulation_configuration(chem_conf['inference_model_configuration']['parameter_list'],
                                                    prior_conf_dict=chem_conf['priors_configuration'],
                                                    photo_ionization_grid=False)

                obj1_model.inference_model()

                # Run the simulation
                # obj1_model.run_sampler(1000, tuning=4000, nchains=6, njobs=6)
                # obj1_model.run_sampler(1000, tuning=1000, nchains=6, njobs=6, init='adapt_diag') # few divergences (75 in one voxel), 23 min
                # obj1_model.run_sampler(1000, tuning=1000, nchains=6, njobs=6, init='advi+adapt_diag' # no divergences 25 min
                # obj1_model.run_sampler(1000, tuning=1000, nchains=6, njobs=6, init='advi+adapt_diag_grad') # no divergences 25 min


                obj1_model.run_sampler(500, tuning=1000, nchains=10, njobs=10, init='advi+adapt_diag') # no divergences 25 min
                obj1_model.save_fit(outputDb, ext_name=chem_ref, output_format='fits')

                # # Load the results
                # fit_pickle = sr.load_MC_fitting(outputDb, ext_name=chem_ref, output_format='fits')
                # inLines = fit_pickle[f'{chem_ref}_inputs'][0]['line_list']
                # inParameters = fit_pickle[f'{chem_ref}_outputs'][0]['parameters_list']
                # inFlux = fit_pickle[f'{chem_ref}_inputs'][0]['line_fluxes']
                # inErr = fit_pickle[f'{chem_ref}_inputs'][0]['line_err']
                # traces_dict = fit_pickle[f'{chem_ref}_traces'][0]

                # # Print the results
                # chemFolder = objFolder/'chemistry'
                # print('-- Model parameters table')
                # figure_file = chemFolder/f'{idx_j}-{idx_i}_table_params'
                # obj1_model.table_mean_outputs(figure_file, inParameters, traces_dict)
                #
                # print('-- Flux values table')
                # figure_file = chemFolder/f'{idx_j}-{idx_i}_table_flux'
                # obj1_model.table_line_fluxes(figure_file, inLines, inFlux, inErr, traces_dict)

                # print('-- Model parameters posterior diagram')
                # figure_file = chemFolder/f'{idx_j}-{idx_i}_plot_params'
                # obj1_model.tracesPosteriorPlot(figure_file, inParameters, traces_dict)
                #
                # print('-- Line flux posteriors')
                # figure_file = chemFolder/f'{idx_j}-{idx_i}_plot_flux'
                # obj1_model.fluxes_distribution(figure_file, inLines, inFlux, inErr, traces_dict)

                # print('-- Model parameters corner diagram')
                # figure_file = chemFolder/f'{idx_j}-{idx_i}_plot_corner'
                # obj1_model.corner_plot(figure_file, inParameters, traces_dict)


    # Show summary
    end = time.time()
    print(f'- Execution time {(end - start)/60:.3f} min')
os.system('spd-say "done"')
