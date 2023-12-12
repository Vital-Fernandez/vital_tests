import numpy as np
import src.specsiser as sr
from delete.data_printing import latex_labels
from pathlib import Path
from astro.data.muse.common_methods import background_color, DARK_PLOT
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS


# Declare data and files location
conf_file = Path('/home/vital/PycharmProjects/vital_tests/astro/data/muse/muse_greenpeas.ini')
seminar_folder = Path('/home/vital/Dropbox/Astrophysics/Seminars/UniVapo 2021/')

obsData = sr.loadConfData(conf_file, group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

voxel_grid_size = obsData['sample_data']['grid_shape_array']

verbose = True

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        objFolder = resultsFolder/obj
        voxelFolder = resultsFolder/obj/'voxel_data'
        db_address = objFolder / f'{obj}_database.fits'
        chem_conf_file = dataFolder/obj/'chemical_model_config.txt'
        fits_results_file = resultsFolder/obj/f'{obj}_chemical.fits'

        # Load chemical conf
        chem_conf = sr.loadConfData(chem_conf_file, group_variables=False)

        # Output files
        chemMaps_fits_address = objFolder/f'{obj}_ChemParamMaps.fits'
        theoLineMaps_fits_address = objFolder/f'{obj}_fitTheoFluxesMaps.fits'

        # ----------------------------------------- Generate the image plots ----------------------------------------
        hdr_plot = fits.getheader(chemMaps_fits_address, extname='PlotConf')
        flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
        halpha_min_level = fits.getval(db_address, keyword=f'P9000', extname=f'H1_6563A_flux')
        halpha_thresd_level = fits.getval(db_address, keyword=f'P8000', extname=f'H1_6563A_flux')

        defaultConf = DARK_PLOT.copy()
        rcParams.update(defaultConf)

        halpha_cmap = cm.gray
        halpha_cmap.set_under(background_color)

        # param lines
        temp_image = fits.getdata(chemMaps_fits_address, 'T_low', ver=1)
        idcs_temp = temp_image > 18000.0
        for i, param in enumerate(['T_low', 'n_e', 'cHbeta', 'NO', 'O']):

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y', 1))

            im = ax.imshow(flux6563_image, cmap=halpha_cmap,
                           norm=colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10))


            ax.set_xlim(100, 210)
            ax.set_ylim(90, 240)

            if param == 'O':
                O2 = np.power(10, fits.getdata(chemMaps_fits_address, 'O2', ver=1)-12)
                O3 = np.power(10, fits.getdata(chemMaps_fits_address, 'O3', ver=1)-12)
                param_image = 12 + np.log10(O2 + O3)
                latex_labels['O'] = r'12+log(O)'
                param_image[idcs_temp] = np.nan
                im2 = ax.imshow(param_image, vmin=7.7, vmax=8.2)

            elif param == 'NO':
                N2 = np.power(10, fits.getdata(chemMaps_fits_address, 'N2', ver=1)-12)
                O2 = np.power(10, fits.getdata(chemMaps_fits_address, 'O2', ver=1)-12)
                param_image = np.log10(N2/O2)
                param_image[idcs_temp] = np.nan
                latex_labels['N/O'] = r'log(N/O)'
                im2 = ax.imshow(param_image)

            else:
                param_image = fits.getdata(chemMaps_fits_address, param, ver=1)
                param_image[idcs_temp] = np.nan
                im2 = ax.imshow(param_image)

            cbar = fig.colorbar(im2, ax=ax)

            param_label = latex_labels[param]
            ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
            plt.tight_layout()
            plt.savefig(seminar_folder/f'{param}_map_region2')
            # plt.show()
