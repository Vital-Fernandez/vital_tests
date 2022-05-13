import numpy as np
import pandas as pd
import src.specsiser as sr
import lime

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import voxel_security_check, muse_grid_sampling_fluxes
from astropy.io import fits

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

# Load the photoionization grid
model_variables = ['logOH', 'logU', 'logNO']
file_address = f'{dataFolder}/HII-CHI-mistry_1Myr_grid_O.txt'

grid_3D = pd.read_csv(file_address, delim_whitespace=True, header=0)
gw = sr.GridWrapper()
grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_3D, axes_columns=model_variables)
grid_interp = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')

# Prepare the fluxes for the sampling
R_v = obsData['Extinction']['R_v']
red_law = obsData['Extinction']['red_law']

# Loop throught the objects and masks
for i, obj in enumerate(objList):

    # Input data
    objFolder = resultsFolder/obj
    db_addresss = objFolder / f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    chemFolder = objFolder / 'chemistry'

    # Output data
    # outputFits = objFolder/f'{obj}_grid_sampling.fits'
    outputFits = objFolder/f'{obj}_grid_sampling_maxErr.fits'

    # Loop throught the line regions
    # for idx_region in [0, 1, 2, 3]:
    for idx_region in [2]:

        # Load fitting configuration
        chem_conf_file = dataFolder / f'grid_sampling_confg_region_{idx_region}.cfg'
        chem_conf = lime.load_cfg(chem_conf_file)
        region_lines = obsData['grid_sampling'][f'region_{idx_region}_line_list']

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
        n_voxels = idcs_voxels.shape[0]
        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            if idx_voxel > 276:

                idx_j, idx_i = idx_pair
                ext_lines = f'{idx_j}-{idx_i}_linelog'
                ext_chem = f'{idx_j}-{idx_i}_gridSampler'
                print(f'\nTreating voxel {idx_j}-{idx_i}: ({idx_voxel}/{n_voxels})')

                # Load voxel fluxes:
                int_series = int_DF.loc[ext_lines]
                err_series = err_DF.loc[ext_lines]

                idcs_obs = ~pd.isnull(int_series) & (int_series.index != 'mask') & (int_series.index.isin(region_lines))
                lineLabels = int_series[idcs_obs].keys().values
                lineInts = int_series[idcs_obs].values
                LineErrs = err_series[idcs_obs].values

                # Define model sampler
                obj1_model = sr.SpectraSynthesizer(grid_sampling=True, grid_interp=grid_interp)
                obj1_model.define_region(lineLabels, lineInts, LineErrs, minErr=np.max(LineErrs/lineInts))
                obj1_model.simulation_configuration(prior_conf_dict=chem_conf['priors_configuration'])
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

    #         # Load voxel lines log
    #         linesLog_BinTable = fits.getdata(fitsLog_address, logLabel, ver=1)
    #         linesDF = Table(linesLog_BinTable).to_pandas()
    #         linesDF.set_index('index', inplace=True)
    #
    #         # Prepare data for HII-CHI-mistry
    #         idcs_inputLines = linesDF.index.isin(labelConver.keys())
    #         input_lines = linesDF.loc[idcs_inputLines].index.values
    #
    #         HII_CHI_mistry_DF = pd.DataFrame()
    #         HII_CHI_mistry_DF.loc[0, 'ID'] = logLabel
    #         flux_Hbeta = linesDF.loc['H1_4861A', 'intg_flux']
    #         for lineLabel in input_lines:
    #             HII_CHI_mistry_DF.loc[0, labelConver[lineLabel]] = linesDF.loc[lineLabel, 'intg_flux'] / flux_Hbeta
    #             HII_CHI_mistry_DF.loc[0, f'e{labelConver[lineLabel]}'] = linesDF.loc[lineLabel, 'intg_err'] / flux_Hbeta
    #         lineSA = HII_CHI_mistry_DF.to_records(index=False) #column_dtypes=default_linelog_types, index_dtypes='<U50')
    #
    #         # Run HII-CHI-mistry
    #         outputSA = epm_HII_CHI_mistry(lineSA, output_file='None', n=200, sed=1, inter=1, HCm_folder=hcm_folder)
    #         linesCol = fits.ColDefs(outputSA)
    #         linesHDU = fits.BinTableHDU.from_columns(linesCol, name=f'{idx_j}-{idx_i}_HIIchimistry')
    #         hdul_lineslog.append(linesHDU)
    #
    #     # Store the drive
    #     hdul_lineslog.writeto(HIIchimistry_fits, overwrite=True, output_verify='fix')
    # end = time.time()
#
# print(f'- Execution time {(end - start)/60:.3f} min')