import numpy as np
import src.specsiser as sr
from src.specsiser.data_printing import latex_labels
from pathlib import Path
from astro.data.muse.common_methods import voxel_security_check, fits_db
from timeit import default_timer as timer
from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS


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
sulfur_bdry = int(obsData['Region_masks']['S3_direct_limit'])
hydrogen_bdry = int(obsData['Region_masks']['H1_direct_limit'])
verbose = True


for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        cube_address = fitsFolder/fileList[i]
        objFolder = resultsFolder/obj
        voxelFolder = resultsFolder/obj/'voxel_data'
        db_addresss = objFolder / f'{obj}_database.fits'
        fitsLog_addresss = objFolder / f'{obj}_linesLog.fits'
        chem_conf_file = dataFolder/obj/'chemical_model_config.txt'
        fits_model_file = resultsFolder/obj/f'{obj}_chemical.fits'

        # Load data
        chem_conf = sr.loadConfData(chem_conf_file, group_variables=False)

        # Load data
        wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')

        # Declare voxels to analyse
        flux6312_image = fits.getdata(db_addresss, f'{ref_flux_line}_flux', ver=1)
        flux6312_levels = np.nanpercentile(flux6312_image, pertil_array)

        flux6563_image = fits.getdata(db_addresss, f'H1_6563A_flux', ver=1)
        flux6563_levels = np.nanpercentile(flux6563_image, pertil_array)

        # Search within that limit
        maFlux_image = np.ma.masked_where((flux6312_image >= flux6312_levels[sulfur_bdry]) &
                                          (flux6563_image > flux6563_levels[hydrogen_bdry]),
                                          flux6563_image)
        idcs_voxels = np.argwhere(maFlux_image.mask)
        n_voxels = idcs_voxels.shape[0]

        if verbose:
            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
            ax.update({'title': r'{} galaxy, $H\alpha$ flux'.format(obj), 'xlabel': r'RA', 'ylabel': r'DEC'})
            im = ax.imshow(maFlux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=flux6563_levels[1],
                           vmin=flux6563_levels[1], base=10))
            cntr1 = ax.contour(flux6312_image, levels=[flux6312_levels[sulfur_bdry]], colors='yellow', alpha=0.5)
            cntr2 = ax.contour(flux6563_image, levels=[flux6563_levels[hydrogen_bdry]], colors='red', alpha=0.5)
            plt.show()

        print(f'\nUsing line [SIII]6312 at percentile {pertil_array[sulfur_bdry]} = {flux6312_levels[sulfur_bdry]:.2f}'
              f' ({idcs_voxels.shape[0]} pixels)')

        parameter_list = ['n_e', 'T_low', 'S3', 'cHbeta']
        parameter_image = np.empty(flux6563_image.shape)

        # Construct the images
        maps_dict = {}
        voxel_n = 0
        with fits.open(fits_model_file) as hdul:

            for parameter in parameter_list:

                parameter_image = np.empty(flux6563_image.shape)
                parameter_image[:] = np.nan

                for idx_voxel, idx_pair in enumerate(idcs_voxels):

                    idx_j, idx_i = idx_pair

                    chem_ref = f'{idx_j}-{idx_i}_chemistry'

                    if chem_ref in hdul:
                        header = hdul[chem_ref].header
                        data = hdul[chem_ref].data
                        trace = data[parameter]
                        parameter_image[idx_j, idx_i] = trace.mean()
                        voxel_n +=1

                maps_dict[parameter] = parameter_image

        # Plot the maps
        print(f'Number of voxels treated: {voxel_n/4}/{n_voxels}')
        for parameter, param_map in maps_dict.items():

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
            im = ax.imshow(maFlux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=flux6563_levels[1],
                                                                              vmin=flux6563_levels[1], base=10))
            im2 = ax.imshow(param_map)
            cbar = fig.colorbar(im2)
            ax.update({'title': r'{} galaxy: {}'.format(obj, latex_labels[parameter]), 'xlabel': r'RA', 'ylabel': r'DEC'})
            plt.tight_layout()
            # plt.savefig(resultsFolder/obj/f'map_{obj}_{parameter}.png', )
            plt.show()
