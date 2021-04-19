import numpy as np
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import voxel_security_check
from timeit import default_timer as timer
from astropy.table import Table
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

ref_flux_line = 'S3_6312A'
int_flux_boundary = -3

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        cube_address = fitsFolder/fileList[i]
        objFolder = resultsFolder/obj
        voxelFolder = resultsFolder/obj/'voxel_data'
        db_addresss = objFolder / f'{obj}_database.fits'
        fitsLog_addresss = objFolder / f'{obj}_linesLog.fits'
        chem_conf_file = dataFolder/obj/'chemical_model_config.txt'

        # Load data
        chem_conf = sr.loadConfData(chem_conf_file, group_variables=False)

        # Declare voxels to analyse
        flux6312_image = fits.getdata(db_addresss, f'{ref_flux_line}_flux', ver=1)
        flux6312_contours = fits.getdata(db_addresss, f'{ref_flux_line}_contour', ver=1)
        flux6312_levels = np.percentile(flux6312_contours[flux6312_contours > 0], pertil_array)
        flux6312_boundary = flux6312_levels[int_flux_boundary]

        flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
        flux6563_contours = fits.getdata(db_addresss, f'H1_6563A_contour', ver=1)
        flux6563_levels = np.percentile(flux6563_contours[flux6563_contours > 0], pertil_array)
        flux6563_boundary = flux6563_levels[int_flux_boundary]

        idcs_bolean = (flux6312_contours >= flux6312_boundary) & (flux6563_contours >= flux6563_boundary)
        idcs_voxels = np.argwhere(idcs_bolean)

        print(f'\nUsing line {ref_flux_line} at perentil {pertil_array[int_flux_boundary]} = {flux6312_boundary:.2f} ({idcs_voxels.shape[0]} pixels)')

        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            idx_j, idx_i = idx_pair

            # Data location
            outputDb = voxelFolder/f'{idx_j}-{idx_i}_fitting.db'
            ext_ref = f'{idx_j}-{idx_i}_linelog'

            # Load voxel data:
            linesDF = Table.read(fitsLog_addresss, ext_ref, character_as_bytes=False).to_pandas()
            linesDF.set_index('index', inplace=True)

            if voxel_security_check:

                # Adjust the lines log for the fitting
                flux_Hbeta = linesDF.loc['H1_4861A', 'intg_flux']

                linesDF.insert(loc=1, column='obsFlux', value=np.nan)
                linesDF.insert(loc=2, column='obsFluxErr', value=np.nan)

                linesDF['obsFlux'] = linesDF['intg_flux']/flux_Hbeta
                linesDF['obsFluxErr'] = linesDF['intg_err']/flux_Hbeta

                # linesDF.rename(columns={'wavelength': 'obsWave'}, inplace=True)
                if ('O2_7319A' in linesDF.index) and ('O2_7330A' in linesDF.index):
                    flux_comb = linesDF.loc['O2_7319A', 'intg_flux'] + linesDF.loc['O2_7330A', 'intg_flux']
                    err_comb = np.sqrt(linesDF.loc['O2_7319A', 'intg_err']**2 + linesDF.loc['O2_7330A', 'intg_err']**2)
                    linesDF.loc['O2_7319A_b'] = None
                    linesDF.loc['O2_7319A_b', ['wavelength', 'obsFlux', 'obsFluxErr']] = 7325, flux_comb/flux_Hbeta, err_comb/flux_Hbeta
                    linesDF.loc['O2_7319A_b', 'ion'] = 'O2'

                # Load emission lines
                input_lines = chem_conf['inference_model_configuration']['input_lines']
                objLinesDF = sr.import_emission_line_data(linesDF=linesDF, include_lines=input_lines)

                # Declare simulation physical properties
                objRed = sr.ExtinctionModel(Rv=chem_conf['simulation_properties']['R_v'],
                                            red_curve=chem_conf['simulation_properties']['reddenig_curve'],
                                            data_folder=chem_conf['data_location']['external_data_folder'])

                objIons = sr.IonEmissivity(tempGrid=chem_conf['simulation_properties']['temp_grid'],
                                           denGrid=chem_conf['simulation_properties']['den_grid'])

                # Generate interpolator from the emissivity grids
                ionDict = objIons.get_ions_dict(np.unique(objLinesDF.ion.values))
                objIons.computeEmissivityGrids(objLinesDF, ionDict, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})

                # Declare chemical model
                objChem = sr.DirectMethod(linesDF=objLinesDF,
                                          highTempIons=chem_conf['simulation_properties']['high_temp_ions_list'])

                # Sampler object
                obj1_model = sr.SpectraSynthesizer()
                obj1_model.define_region(objLinesDF, objIons, objRed, objChem)

                obj1_model.simulation_configuration(chem_conf['inference_model_configuration']['parameter_list'],
                                                    prior_conf_dict=chem_conf['priors_configuration'],
                                                    photo_ionization_grid=False)
                obj1_model.inference_model()

                # Run the simulation
                obj1_model.run_sampler(outputDb, 3000, tuning=3000, nchains=3, njobs=3)

                # Plot the results
                fit_results = sr.load_MC_fitting(outputDb)

                print('-- Printing results')
                figure_file = voxelFolder / f'{idx_j}-{idx_i}_MeanOutputs'
                obj1_model.table_mean_outputs(figure_file, fit_results)

                figure_file = voxelFolder / f'{idx_j}-{idx_i}_FluxComparison'
                obj1_model.table_line_fluxes(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})

                figure_file = voxelFolder / f'{idx_j}-{idx_i}_ParamsPosteriors.png'
                obj1_model.tracesPosteriorPlot(figure_file, fit_results)

                figure_file = voxelFolder / f'{idx_j}-{idx_i}_lineFluxPosteriors.png'
                obj1_model.fluxes_distribution(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})

                figure_file = voxelFolder / f'{idx_j}-{idx_i}_cornerPlot.png'
                obj1_model.corner_plot(figure_file, fit_results)
                obj1_model.savefig(figure_file, resolution=200)

            #     # Declare sampling properties
            #     obj1_model.simulation_configuration(chem_conf['inference_model_configuration']['parameter_list'],
            #                                         prior_conf_dict=chem_conf['priors_configuration'])
            #
            #     # Declare simulation inference model
            #     obj1_model.inference_model(include_Thigh_prior=chem_conf['inference_model_configuration']['T_high_check'])
            #
            #     # Run the simulation
            #     obj1_model.run_sampler(outputDb, 5000, 2000, njobs=1)
            #
            #     if outputDb.is_file() and Path(outputDb).with_suffix('.txt').is_file():
            #
            #         # Plot the results
            #         fit_results = sr.load_MC_fitting(outputDb)
            #
            #         # Print the results
            #         print('-- Model parameters table')
            #         figure_file = voxelFolder / f'{idx_j}-{idx_i}_MeanOutputs'
            #         obj1_model.table_mean_outputs(figure_file, fit_results)
            #
            #         print('-- Flux values table')
            #         figure_file = voxelFolder / f'{idx_j}-{idx_i}_FluxComparison'
            #         obj1_model.table_line_fluxes(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})
            #
            #         print('-- Model parameters posterior diagram')
            #         figure_file = voxelFolder / f'{idx_j}-{idx_i}_ParamsPosteriors.png'
            #         obj1_model.tracesPosteriorPlot(figure_file, fit_results)
            #
            #         print('-- Line flux posteriors')
            #         figure_file = voxelFolder / f'{idx_j}-{idx_i}_lineFluxPosteriors.png'
            #         obj1_model.fluxes_distribution(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})

