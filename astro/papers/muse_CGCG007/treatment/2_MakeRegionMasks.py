import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import lineAreas, store_frame_to_fits
from astropy.io import fits
from matplotlib import pyplot as plt, cm, colors
from astropy.wcs import WCS

# Declare data and files location
obsData = sr.loadConfData('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

dict_errs = {}
dict_nan_values = {}

ref_flux_line = 'S3_6312A'

for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_addresss = objFolder/f'{obj}_database.fits'
    mask_address = dataFolder/f'{obj}_mask.txt'

    # Output data:
    fitsLog_addresss = objFolder/f'{obj}_linesLog.fits'
    hdul_lineslog = fits.HDUList()

    # Load data
    wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')
    wave_rest = wave / (1 + z_objs[i])
    mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)

    # Declare voxels to analyse
    invers_pertil_array = pertil_array[::-1]
    flux6312_image = fits.getdata(db_addresss, f'{ref_flux_line}_flux', ver=1)
    flux6312_levels = np.nanpercentile(flux6312_image, invers_pertil_array)

    flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
    flux6563_levels = np.nanpercentile(flux6563_image, invers_pertil_array)

    Halpha_idxMinPercentil = 6
    Halpha_min_level = flux6563_levels[5]

    SIII_idxMinPercentil = 3

    for idx_contour in np.arange(Halpha_idxMinPercentil):

        # Search within that limit
        if idx_contour < SIII_idxMinPercentil:

            if idx_contour == 0:
                maFlux_image = np.ma.masked_where((flux6312_image >= flux6312_levels[idx_contour]) &
                                                  (flux6563_image > Halpha_min_level),
                                                  flux6563_image)
            else:
                maFlux_image = np.ma.masked_where((flux6312_image >= flux6312_levels[idx_contour]) &
                                                  (flux6312_image < flux6312_levels[idx_contour - 1]) &
                                                  (flux6563_image > Halpha_min_level),
                                                  flux6563_image)
        elif idx_contour == 3:
            maFlux_image = np.ma.masked_where((flux6563_image > flux6563_levels[idx_contour]) &
                                              (flux6312_image < flux6312_levels[idx_contour - 1]) &
                                              (flux6563_image > Halpha_min_level),
                                              flux6563_image)

        else:
            maFlux_image = np.ma.masked_where((flux6563_image >= flux6563_levels[idx_contour]) &
                                              (flux6563_image < flux6563_levels[idx_contour-1]) &
                                              (flux6563_image > Halpha_min_level),
                                              flux6563_image)


        idcs_voxels = np.argwhere(maFlux_image.mask)
        print(f'region {idx_contour} ({idcs_voxels.shape[0]} pixels)')

        # if verbose:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
        ax.update({'title': r'{} galaxy, $H\alpha$ flux'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
        im = ax.imshow(maFlux_image, cmap=cm.gray,
                       norm=colors.SymLogNorm(linthresh=flux6563_levels[-2], vmin=flux6563_levels[-2], base=10))
        # plt.savefig(objFolder/f'region_{idx_contour}_mask')
        plt.show()

        mask_name = f'region_{idx_contour}'
        mask_hdu = fits.ImageHDU(name=mask_name, data=maFlux_image.mask.astype(int), ver=1)
        store_frame_to_fits(db_addresss, fits_hdu=mask_hdu, ext_name=mask_name)

