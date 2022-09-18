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
# ref_simulations = ['localErr', 'HIICHImistry', 'maxErr',  'minOneErr']
ref_simulations = ['HIICHImistry', 'maxErr']
i_start = {'minOneErr': 264}

tech_label = 'GridSampling'

# Load the photoionization grid
model_variables = ['logOH', 'logU', 'logNO']

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
        # for idx_region in [0, 1, 2]:
        for idx_region in [0, 1]:

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
                obj1_model.save_fit(outputFits, ext_chem, output_format='fits')
