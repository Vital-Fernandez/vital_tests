import numpy as np
import src.specsiser as sr
import time
import lime

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import voxel_security_check, chemical_lines_indexing
from astropy.io import fits

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
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

R_v = obsData['Extinction']['R_v']
red_law = obsData['Extinction']['red_law']

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = f'{objFolder}/chemistry'
    obsLog_addresss = objFolder / f'{obj}_linesLog.fits'
    absLog_addresss = objFolder / f'{obj}_linesLog_abs.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Output files
    outputDb = objFolder/f'{obj}_chemical.fits'

    # Loop throught the line regions
    start = time.time()
    # for idx_region in [0, 1, 2, 3]:
    for idx_region in [2]:

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
        input_lines = chem_conf['inference_model_configuration']['input_lines']
        emis_grid_interp = sr.emissivity_grid_calc(lines_array=input_lines, comp_dict=merge_dict)

        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            idx_j, idx_i = idx_pair

            print(f'\nTreating voxel {idx_j}-{idx_i}: ({idx_voxel}/{n_voxels})')

            # Data location
            chem_ref = f'{idx_j}-{idx_i}_chemistry'

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

            if voxel_security_check(linesDF):

                # Declare object
                obj1_model = sr.SpectraSynthesizer(emis_grid_interp)

                # Declare simulation physical properties
                obj1_model.define_region(lineLabels, lineFluxes, lineErr, lineflambda)

                # Sampling configuration
                obj1_model.simulation_configuration(prior_conf_dict=chem_conf['priors_configuration'],
                                                    highTempIons=chem_conf['simulation_properties']['high_temp_ions_list'])
                # Theoretical model
                obj1_model.inference_model()

                obj1_model.run_sampler(500, 2000, nchains=10, njobs=10, init='advi+adapt_diag')
                obj1_model.save_fit(outputDb, ext_name=chem_ref, output_format='fits')

            #     # Load the results
            #     print(chem_ref)
            #     fit_pickle = sr.load_fit_results(outputDb, ext_name=chem_ref, output_format='fits')
            #     inLines = fit_pickle[f'{chem_ref}_inputs'][0]['line_list']
            #     inParameters = fit_pickle[f'{chem_ref}_outputs'][0]['parameters_list']
            #     inFlux = fit_pickle[f'{chem_ref}_inputs'][0]['line_fluxes']
            #     inErr = fit_pickle[f'{chem_ref}_inputs'][0]['line_err']
            #     traces_dict = fit_pickle[f'{chem_ref}_traces'][0]
            #
            #     # Print the results
            #     print('-- Model parameters table')
            #     figure_file = f'{chemFolder}/{chem_ref}_fitted_fluxes'
            #     sr.table_fluxes(figure_file, inLines, inFlux, inErr, traces_dict, merge_dict)

            # # Print the results
            # print('-- Fitted fluxes table')
            # figure_file = f'{chemFolder}/{chem_ref}_MeanOutputs'
            # sr.table_params(figure_file, inParameters, traces_dict)
            #
            # print('-- Model parameters posterior diagram')
            # figure_file = f'{chemFolder}/{chem_ref}_traces_plot.png'
            # sr.plot_traces(figure_file, inParameters, traces_dict)
            #
            # print('-- Line flux posteriors')
            # figure_file = f'{chemFolder}/{chem_ref}_fluxes_grid.png'
            # sr.plot_flux_grid(figure_file, inLines, inFlux, inErr, traces_dict)
            #
            # print('-- Model parameters corner diagram')
            # figure_file = f'{chemFolder}/{chem_ref}_corner.png'
            # sr.plot_corner(figure_file, inParameters, traces_dict)


    # # Show summary
    # end = time.time()
    # print(f'- Execution time {(end - start)/60:.3f} min')
