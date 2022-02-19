import numpy as np
import pandas as pd
import src.specsiser as sr
import time
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

# Loop throught the objects and masks
for i, obj in enumerate(objList):

    # Input data
    objFolder = resultsFolder/obj
    db_addresss = objFolder / f'{obj}_database.fits'
    fitsLog_addresss = objFolder / f'{obj}_linesLog.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    outputDb = objFolder/f'{obj}_chemical.fits'
    chemFolder = objFolder / 'chemistry'

    # Output data
    outputFits = objFolder/f'{obj}_grid_fit.fits'

    # Loop throught the line regions
    start = time.time()
    for idx_region in [0]:
    # for idx_region in [0, 1, 2]:

        # Voxel mask
        region_label = f'mask_{idx_region}'
        region_mask = fits.getdata(maskFits_address, region_label, ver=1)
        region_mask = region_mask.astype(bool)
        idcs_voxels = np.argwhere(region_mask)

        # Simulation time estimation
        n_voxels = idcs_voxels.shape[0]
        print(f'\nTreating {region_label} consisting of {n_voxels} voxels')

        # Region chemical configuration
        chem_conf_file = dataFolder / f'HII_CHIM_TRY_config_region_{idx_region}.cfg'
        chem_conf = lime.load_cfg(chem_conf_file)

        # Load emission lines
        input_lines = chem_conf['inference_model_configuration']['input_lines_list']
        ion_array, wave_array, latex_array = lime.label_decomposition(input_lines)

        # Loop through the region voxels
        n_voxels = idcs_voxels.shape[0]
        for idx_voxel, idx_pair in enumerate(idcs_voxels):

            idx_j, idx_i = idx_pair
            print(f'\nTreating voxel {idx_j}-{idx_i}: ({idx_voxel}/{n_voxels})')

            # Load voxel fluxes:
            ext_lines = f'{idx_j}-{idx_i}_linelog'
            linesDF = lime.load_lines_log(fitsLog_addresss, ext_lines)

            # Load the results
            ext_chem = f'{idx_j}-{idx_i}_chemistry'
            fit_results = sr.load_fit_results(outputDb, ext_name=ext_chem, output_format='fits')
            cHBeta = fit_results[f'{ext_chem}_outputs'][1]['cHbeta']
            print(f'- cHBeta: {cHBeta}')

            # Prepare the fluxes for the sampling
            Rv = chem_conf['simulation_properties']['R_v']
            redLaw = chem_conf['simulation_properties']['reddenig_curve']
            input_log = muse_grid_sampling_fluxes(input_lines, linesDF, cHBeta, Rv, redLaw)

            # Define model sampler
            ext_chem = f'{idx_j}-{idx_i}_gridSampler'
            obj1_model = sr.SpectraSynthesizer(grid_sampling=True, grid_interp=grid_interp)
            obj1_model.define_region(input_log.index.values, input_log.obsInt.values, input_log.obsIntErr.values, lineFlambda=np.zeros(len(input_log)))
            obj1_model.simulation_configuration(prior_conf_dict=chem_conf['priors_configuration'])
            obj1_model.photoionization_sampling(model_variables)
            obj1_model.run_sampler(1000, 2000, nchains=8, njobs=8, init='advi')
            # obj1_model.run_sampler(1000, 2000, nchains=8, njobs=8, init='advi')
            obj1_model.save_fit(outputFits, ext_chem, output_format='fits')

            # Load the results
            fit_results = sr.load_fit_results(outputFits, ext_name=ext_chem, output_format='fits')
            inLines = fit_results[f'{ext_chem}_inputs'][0]['line_list']
            inParameters = fit_results[f'{ext_chem}_outputs'][0]['parameters_list']
            inFlux = fit_results[f'{ext_chem}_inputs'][0]['line_fluxes']
            inErr = fit_results[f'{ext_chem}_inputs'][0]['line_err']
            traces_dict = fit_results[f'{ext_chem}_traces'][0]

            print('-- Model parameters posterior diagram')
            figure_file = f'{chemFolder}/{ext_chem}_trace_plot.png'
            sr.plot_traces(figure_file, inParameters, traces_dict)



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