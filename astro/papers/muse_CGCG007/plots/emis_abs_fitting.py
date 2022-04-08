import numpy as np
import pandas as pd
import time
import lime
from scipy.interpolate import interp1d
from pathlib import Path
from astropy.io import fits

from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits
from progressbar import progressbar

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
thres_array = obsData['sample_data']['detect_lim_array']

dict_errs = {}
dict_nan_values = {}

verbose = True

for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_addresss = objFolder/f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    mask_address = dataFolder / f'{obj}_region1_mask.txt'

    # Load data
    wave, cube, header = import_muse_fits(cube_address)
    mask = lime.load_lines_log(mask_address)

    idx_j, idx_i = (177, 149)

    flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
    flux_err = np.sqrt(cube[:, idx_j, idx_i].var.data) * norm_flux

    voxel = lime.Spectrum(wave, flux_voxel, input_err=flux_err, redshift=z_objs[i], norm_flux=norm_flux)
    # voxel.plot_spectrum()

    line = 'H1_4861A'
    fit_conf = obsData['emis_abs_line_fitting']
    mask.loc[line, 'w3'] = mask.loc[line, 'w2']
    mask.loc[line, 'w4'] = mask.loc[line, 'w5']
    mask.loc[line, 'w2'] = mask.loc[line, 'w2'] - 5
    mask.loc[line, 'w5'] = mask.loc[line, 'w5'] + 5
    mask_array = mask.loc[line, 'w1':'w6'].values

    voxel.fit_from_wavelengths('H1_4861A_b', mask_array, fit_conf, emission=True)
    voxel.display_results(fit_report=True)

    # line = 'H1_4861A_b'
    #
    # voxel.fit_from_wavelengths(line, mask_array, fit_conf, fit_method='least_squares')
    # voxel.display_results(fit_report=True)