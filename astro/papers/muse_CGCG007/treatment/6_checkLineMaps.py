import numpy as np
import lime
from pathlib import Path
from astropy.io import fits
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits, param_images

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
norm_flux = obsData['sample_data']['norm_flux']

for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_address = objFolder / f'{obj}_database.fits'
    fitsLog_address = objFolder / f'{obj}_linesLog.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Labels for the axes
    wave, cube, header = import_muse_fits(cube_address)
    flux = cube.data.data * norm_flux

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
    print(fits.info(fitsLog_address))
    # Load data
    ax_conf = {'image': {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'MUSE CGCG007-025'}}
    lime.CubeFitsInspector(wave, flux, Halpha_image, SIII_image, SIII_contourLevels, fits_header=header,
                           lines_log_address=fitsLog_address, axes_conf=ax_conf, redshift=z_objs[i])

