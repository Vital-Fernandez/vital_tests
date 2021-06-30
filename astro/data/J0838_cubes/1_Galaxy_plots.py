import numpy as np
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors
from src.specsiser.print.plot import STANDARD_PLOT
from astropy.io import fits
from src.specsiser.tools.line_fitting import EmissionFitting


# Declare data and files location
obsConf = sr.loadConfData('J0838_cubes.ini')
fitsFolder = Path(obsConf['data_location']['fits_folder'])
dataFolder = Path(obsConf['data_location']['data_folder'])
resultsFolder = Path(obsConf['data_location']['results_folder'])

fileList = obsConf['data_location']['file_list']
objList = obsConf['data_location']['object_list']
z_list = obsConf['sample_data']['z_array']
norm_flux = obsConf['sample_data']['norm_flux']
percentil_array = obsConf['sample_data']['percentil_array']

# Plot set up
labelsDict = {'xlabel': r'RA',
              'ylabel': r'DEC'}
defaultConf = STANDARD_PLOT.copy()
defaultConf.update(labelsDict)
rcParams.update({})

ef = EmissionFitting()

lines_cube = {objList[0]: ['O3_5007A'],
              objList[1]: ['H1_6563A']}

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        objFolder = resultsFolder
        cube_address_i = fitsFolder/fileList[i]
        mask_address = dataFolder/obsConf['data_location']['mask_global']

        # Output data
        db_address = objFolder / f'{obj}_database.fits'

        # Load the data
        wave, data, header = sr.import_fits_data(cube_address_i, instrument='fits-cube', frame_idx=0)
        mask_global_DF = sr.lineslogFile_to_DF(mask_address)

        # Create empty fits file
        hdul_log = fits.HDUList()
        hdul_log.append(fits.PrimaryHDU())

        # Second page for the fits file plot configuration
        col_waves = fits.Column(name='wave', array=wave, format='1E')
        hdu_table = fits.BinTableHDU.from_columns([col_waves], name='PlotConf', header=header)
        hdul_log.append(hdu_table)
        # hdul_log[1].header = header

        # Create flux maps for the main lines:
        for lineLabel in lines_cube[obj]:

            lineWaves = mask_global_DF.loc[lineLabel, 'w1':'w6']
            ion, wavelength, latexLabel = sr.label_decomposition(lineLabel, scalar_output=True)

            idcsLineRegion, idcsContRegion = ef.define_masks(wave, data, lineWaves)
            line_slice, line_cont = data[idcsLineRegion, :, :], data[idcsContRegion, :, :]

            avg = np.mean(line_slice, axis=0)
            rms = np.sqrt(np.mean(np.power(line_slice - avg, 2), axis=0))
            SNR_image = avg / rms

            flux_image = line_slice.sum(axis=0)
            levelContours = np.nanpercentile(flux_image, percentil_array)

            plot_image_file = objFolder/f'{obj}_{lineLabel}_contours.png'

            # # Extract cube slice using mpdaf defult tools.
            # # This requires the input wavelengths to be on the same scale as in the cube
            # line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
            # flux_image = line_image.data.data
            # levelContours = np.nanpercentile(flux_image, pertil_array)

            # Store fluxes and contours
            hdu_image = fits.ImageHDU(name=f'{lineLabel}_flux', data=flux_image, ver=1)
            for idx, level in enumerate(levelContours):
                level_label = f'hierarch P{int(percentil_array[idx]*100)}'
                print(level_label)
                hdu_image.header[level_label] = level
            hdul_log.append(hdu_image)

            # Plot the image:
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot()

            im = ax.imshow(flux_image, cmap=cm.gray)
            cntr1 = ax.contour(flux_image, levels=levelContours, cmap='viridis')


            # im = ax.imshow(flux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=levelContours[-3], vmin=levelContours[-3], base=10))
            # cntr1 = ax.contour(flux_image, levels=levelContours[::-1][idx_ion_boundary:], cmap='viridis', norm=colors.LogNorm())
            #
            # for idx, percentile in enumerate(pertil_array[::-1][idx_ion_boundary:]):
            #     label = r'$P_{{{}}}$({})'.format(idx, latexLabel)
            #     cntr1.collections[idx].set_label(label)
            # ax.legend()
            #
            ax.update({'title': r'{} galaxy, {} flux'.format(obj, latexLabel), 'xlabel': r'Voxel', 'ylabel': r'Voxel'})
            plt.savefig(plot_image_file, bbox_inches='tight')
            plt.show()

        # Store the drive
        hdul_log.writeto(db_address, overwrite=True, output_verify='fix')

    # hdr = fits.getheader(db_addresss, extname='PlotConf')
    # for lineLabel, lineLimits in lineAreas.items():
    #     flux_image = fits.getdata(db_addresss, f'{lineLabel}_flux', ver=1)
    #     flux_contours = fits.getdata(db_addresss, f'{lineLabel}_contour', ver=1)
    #
    #     # Define image countours based on the flux percentiles
    #     levelFlux_i = np.percentile(flux_contours[flux_contours > 0], pertil_array)
    #     levels_text_i = ['None'] * len(levelFlux_i)
    #     for idx, per in enumerate(pertil_array):
    #         levels_text_i[idx] = f'{levelFlux_i[idx]:.2f} $P_{{{per}}}$ $log(F_{{\lambda}})$'
    #
    #     # Plot the image:
    #     fig = plt.figure(figsize=(12, 8))
    #     ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y', 1))
    #
    #     frame_size = flux_image.shape
    #     x, y = np.arange(0, frame_size[1]), np.arange(0, frame_size[0])
    #     X, Y = np.meshgrid(x, y)
    #
    #     CS3 = ax.contourf(X, Y, flux_contours, levels=levelFlux_i)
    #     cbar = fig.colorbar(CS3)
    #     cbar.ax.set_yticklabels(levels_text_i)
    #     ax.set_facecolor('black')
    #     ax.update(labelsDict)
    #     ax.set_title(f'Galaxy {obj} {lineLabel}')
    #     plt.show()
