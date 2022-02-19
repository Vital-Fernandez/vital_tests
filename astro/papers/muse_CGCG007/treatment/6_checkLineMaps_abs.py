import numpy as np
import lime
from pathlib import Path
from astropy.io import fits
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_fado_cube, import_muse_fits

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList_obs = obsData['data_location']['file_list']
fileList_abs = obsData['data_location']['file_abs_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
fitsFolder_abs = Path(obsData['data_location']['fits_abs_folder'])

z_objs = obsData['sample_data']['z_array']
norm_flux_abs = obsData['sample_data']['norm_flux_abs']
norm_flux = obsData['sample_data']['norm_flux']

for i, obj in enumerate(objList):

    # Data location
    cube_address_obs = fitsFolder/fileList_obs[i]
    cube_address_abs = fitsFolder_abs/fileList_abs[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_address = objFolder / f'{obj}_database.fits'
    fitsLog_address = objFolder / f'{obj}_linesLog_abs.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Open emission cube
    wave, cube, header = import_muse_fits(cube_address_obs)
    flux = cube.data.data * norm_flux

    # Open absorption cube
    wave, cube_abs, header_abs = import_fado_cube(cube_address_abs)

    # Establish the line intervals for the cube flux image slices
    line_regions = {'H1_6563A': np.array([6528, 6591]) * (1 + z_objs[i]),
                    'S3_6312A': np.array([6308, 6317]) * (1 + z_objs[i])}

    # Use Halpha for the background
    idcs_Halpha = np.searchsorted(wave, line_regions['H1_6563A'])
    Halpha_image = flux[idcs_Halpha[0]:idcs_Halpha[1], :, :].sum(axis=0)

    # Use SII lines for the foreground contours
    idcs_line = np.searchsorted(wave, line_regions['S3_6312A'])
    SIII_image = flux[idcs_line[0]:idcs_line[1], :, :].sum(axis=0)

    # Establishing the countours by percentiles
    percentil_array = np.array([99, 99.9])
    SIII_contourLevels = np.nanpercentile(SIII_image, percentil_array)

    # Load data
    ax_conf = {'image': {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'MUSE CGCG007-025'}}
    lime.CubeFitsInspector(wave, cube_abs, Halpha_image, SIII_image, SIII_contourLevels, fits_header=header_abs,
                           lines_log_address=fitsLog_address, axes_conf=ax_conf)

