import numpy as np
import pandas as pd
import time
import lime

from pathlib import Path
from astropy.io import fits

pd.set_option('display.max_columns', None)

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

R_v = obsData['Extinction']['R_v']
red_law = obsData['Extinction']['red_law']

HII_chim_folder = Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v5.22')

image_size = obsData['sample_data']['grid_shape_array'].astype(int)

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = objFolder/'chemistry'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    db_address = objFolder / f'{obj}_database.fits'

    # Header data with the cosmological coordinates
    hdr_plot = fits.getheader(db_address, extname='PlotConf')

    plot_dict = {'CRPIX1': hdr_plot['CRPIX1'],
                 'CRPIX2': hdr_plot['CRPIX2'],
                 'CD1_1': hdr_plot['CD1_1'],
                 'CD1_2': hdr_plot['CD1_2'],
                 'CD2_1': hdr_plot['CD2_1'],
                 'CD2_2': hdr_plot['CD2_2'],
                 'CUNIT1': hdr_plot['CUNIT1'],
                 'CUNIT2': hdr_plot['CUNIT2'],
                 'CTYPE1': hdr_plot['CTYPE1'],
                 'CTYPE2': hdr_plot['CTYPE2']}
    ext_hdr = fits.Header(plot_dict)

    # Get mask indeces:
    mask_list = ['MASK_0', 'MASK_1', 'MASK_2']
    spatial_mask_dict = {}
    with fits.open(maskFits_address) as hdu_masks:
        for mask_name in mask_list:
            mask_data = hdu_masks[mask_name].data.astype(bool)
            spatial_mask_dict[mask_name] = mask_data
    total_mask = np.array(list(spatial_mask_dict.values()))
    total_mask = total_mask.sum(axis=0).astype(bool)

    # Storing dictionary
    param_list = ['OH', 'logU', 'NO']
    param_maps_dict = dict.fromkeys(param_list + ['OH_err', 'logU_err', 'NO_err'], None)
    for item in param_maps_dict:
        param_maps_dict[item] = np.full(image_size, np.nan)

    # Loop throught the line regions
    # for idx_region in [0, 1, 2, 3]:
    for idx_region in [0, 1, 2]:

        fit_output_file = HII_chim_folder/f'region_{idx_region}_intensities_hcm-output.dat'

        if fit_output_file.is_file():

            # Open HII-CHI-mistry file
            print(f'- Converting region {idx_region}')
            output_DF = pd.read_csv(fit_output_file, delim_whitespace=True, skiprows=7, index_col=0)

            # Loop through the mask voxels
            idcs_voxels = np.argwhere(spatial_mask_dict[f'MASK_{idx_region}'])
            for idx_pair in idcs_voxels:
                idx_j, idx_i = idx_pair
                cord_label = f'{idx_j}-{idx_i}_line'
                for param in param_list:
                    param_value = output_DF.loc[cord_label, param]
                    param_err = output_DF.loc[cord_label, f'e{param}']
                    print(param, param_value, param_err)
                    param_maps_dict[param][idx_j, idx_i] = param_value
                    param_maps_dict[f'{param}_err'][idx_j, idx_i] = param_err

    # Save as fits files
    for param in param_list:

        # Fits extensions
        idcs_data = ~np.isnan(param_maps_dict[param])

        paramHDUs = fits.HDUList()
        paramHDUs.append(fits.PrimaryHDU())
        paramHDUs.append(fits.ImageHDU(name=param, data=param_maps_dict[param], header=ext_hdr, ver=1))
        paramHDUs.append(fits.ImageHDU(name=f'{param}_err', data=param_maps_dict[f'{param}_err'], header=ext_hdr, ver=1))

        # Extension files
        output_file = chemFolder/f'{param}_HII-CHI-mistry_output.fits'
        paramHDUs.writeto(output_file, overwrite=True, output_verify='fix')

