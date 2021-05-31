import numpy as np
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import voxel_security_check, fits_db
from timeit import default_timer as timer
from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
import time


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

verbose = True

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        objFolder = resultsFolder/obj
        voxelFolder = resultsFolder/obj/'voxel_data'
        db_addresss = objFolder / f'{obj}_database.fits'
        fitsLog_addresss = objFolder / f'{obj}_linesLog.fits'

        # Output files
        outputFitsFile = resultsFolder/obj/f'{obj}_chemical.fits'

        # Loop throught the line regions
        start = time.time()
        for idx_region in [0, 1, 2]:

            # Voxel mask
            region_label = f'region_{idx_region}'
            region_mask = fits.getdata(db_addresss, region_label, ver=1)
            region_mask = region_mask.astype(bool)
            idcs_voxels = np.argwhere(region_mask)
            n_voxels = idcs_voxels.shape[0]
            print(f'\nTreating {region_label} consisting of {n_voxels}. The estimated time is {(n_voxels*1.85)/60:0.1f} hrs')

            # Region chemical configuration
            chem_conf_file = dataFolder/obj/f'chemical_model_region_{idx_region}.txt'
            chem_conf = sr.loadConfData(chem_conf_file, group_variables=False)

            for idx_voxel, idx_pair in enumerate(idcs_voxels):

                idx_j, idx_i = idx_pair
                print(f'\nTreating voxel {idx_j}-{idx_i}: ({idx_voxel}/{n_voxels})')

                # Data location
                outputDb = voxelFolder/f'{idx_j}-{idx_i}_fitting.db'
                chem_ref = f'{idx_j}-{idx_i}_chemistry'

                # Load voxel data:
                ext_ref = f'{idx_j}-{idx_i}_linelog'
                linesDF = Table.read(fitsLog_addresss, ext_ref, character_as_bytes=False).to_pandas()
                linesDF.set_index('index', inplace=True)
                if voxel_security_check:

                    linesDF.insert(loc=1, column='obsFlux', value=np.nan)
                    linesDF.insert(loc=2, column='obsFluxErr', value=np.nan)

                    idcs_gaussian = (linesDF.blended != 'None') & (~linesDF.index.str.contains('_b'))
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

                    # Load emission lines
                    input_lines = chem_conf['inference_model_configuration']['input_lines']
                    objLinesDF = sr.import_emission_line_data(linesDF=linesDF, include_lines=input_lines)

                    if 'S3_6312A' in objLinesDF.index:

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
                        obj1_model.run_sampler(outputDb, 2000, tuning=3000, nchains=4, njobs=4)

                        # Plot the results
                        fit_results = sr.load_MC_fitting(outputDb)

                        # Store the results
                        fits_db(outputFitsFile, fit_results, chem_ref)

                        # print('-- Printing results')
                        # figure_file = voxelFolder / f'{idx_j}-{idx_i}_MeanOutputs'
                        # obj1_model.table_mean_outputs(figure_file, fit_results)
                        #
                        # figure_file = voxelFolder / f'{idx_j}-{idx_i}_FluxComparison'
                        # obj1_model.table_line_fluxes(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})
                        #
                        # figure_file = voxelFolder / f'{idx_j}-{idx_i}_ParamsPosteriors.png'
                        # obj1_model.tracesPosteriorPlot(figure_file, fit_results)
                        #
                        # figure_file = voxelFolder / f'{idx_j}-{idx_i}_lineFluxPosteriors.png'
                        # obj1_model.fluxes_distribution(figure_file, fit_results, combined_dict={'O2_7319A_b': 'O2_7319A-O2_7330A'})
                        #
                        # figure_file = voxelFolder / f'{idx_j}-{idx_i}_cornerPlot.png'
                        # obj1_model.corner_plot(figure_file, fit_results)

        # Show summary
        end = time.time()
        print(f'- Execution time {(end - start)/60:.3f} min')
