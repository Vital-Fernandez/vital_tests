import numpy as np
import src.specsiser as sr
from pathlib import Path
from mpdaf.obj import Cube
from matplotlib import pyplot as plt, rcParams
import astropy.units as u
from mpdaf.obj import deg2sexa
from astropy.wcs import WCS
from src.specsiser.print.plot import STANDARD_PLOT
from astropy.visualization import mpl_normalize, SqrtStretch

# Declare data and files location
conf_file_address = '../muse_greenpeas.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)
objList = obsData['sample_data']['object_list']
fileList = obsData['sample_data']['file_list']
dataFolder = obsData['sample_data']['data_folder']
idx_voxel = (170, 170)
lineAreas = {'H1_6563A': (6558.0, 6568.0), 'S3_6312A': (6308.15, 6317.25), 'O3_5007A': (5002.0, 5013.0), 'S2_6717A': (6717.0, 6734.0)}
z_objs = obsData['sample_data']['z_array']


for i, obj in enumerate(objList):

        # Load the data
        print(f'\n- {obj}')
        file_address_i = f'{dataFolder}/{fileList[i]}'
        wave, cube, header = sr.import_fits_data(file_address_i, instrument='MUSE')
        wave = wave / (1 + z_objs[i])
        print(f'\n-- {header["OBJECT"]}')

        # Get astronomical coordinates one pixel
        coord_sky = cube.wcs.pix2sky(idx_voxel, unit=u.deg)
        dec, ra = deg2sexa(coord_sky)[0]
        wcs_cube = WCS(cube.data_header)

        # Treat all the lines
        for lineLabel, lineLimits in lineAreas.items():

            # Get line flux region
            ion, lineWave, latexlabel = sr.label_decomposition([lineLabel])
            # lineIdcs = np.searchsorted(wave, np.array(lineLimits))
            # lineImage = cube[lineIdcs[0]:lineIdcs[1], :, :].sum(axis=0)/wave[lineIdcs[0]:lineIdcs[1]].size
            # print(lineLabel, wave[lineIdcs[0]:lineIdcs[1]].size)
            # flux_image = lineImage.data.data

            line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
            flux_image = line_image.data.data

            # Plot line image map with coordinates
            labelsDict = {'xlabel': r'RA',
                          'ylabel': r'DEC',
                          'title': r'Galaxy {} {}'.format(obj, latexlabel[0])}

            # Plot Configuration
            defaultConf = STANDARD_PLOT.copy()
            defaultConf.update(labelsDict)
            rcParams.update({})

            # Selecting plotting value pixels
            frame_size = flux_image.shape
            x, y = np.arange(0, frame_size[1]), np.arange(0, frame_size[0])
            X, Y = np.meshgrid(x, y)

            log_flux = np.log10(flux_image)
            flux_contours = np.zeros(flux_image.shape)
            idcs_removed = np.logical_or(log_flux < 0.0, np.isnan(log_flux))
            flux_contours[~idcs_removed] = log_flux[~idcs_removed]

            percent_list = np.array([0, 80, 90, 95, 99.5, 99.90, 99.99])
            levels = np.percentile(flux_contours, percent_list)
            levels_text = ['None'] * levels.size
            for idx, per in enumerate(percent_list):
                levels_text[idx] = f'{levels[idx]:.2f} $P_{{{per}}}$'

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot(projection=wcs_cube, slices=('x', 'y', 1))

            CS3 = ax.contourf(X, Y, flux_contours, levels=levels)
            cbar = fig.colorbar(CS3)
            cbar.ax.set_yticklabels(levels_text)

            # CS3 = ax.contourf(X, Y, flux_contours)
            # cbar = fig.colorbar(CS3)

            # ax.set_facecolor('black')
            ax.update(labelsDict)
            imageName = fileList[i].replace('.fits', f'_{lineLabel}_contours.png')
            mapAddress = file_address_i = f'{dataFolder}/{imageName}'
            # plt.savefig(mapAddress, bbox_inches='tight')
            plt.show()







# import numpy as np
# import src.specsiser as sr
# from pathlib import Path
# from mpdaf.obj import Cube
# import matplotlib.pyplot as plt
# import astropy.units as u
# from mpdaf.obj import deg2sexa
# from astropy.wcs import WCS
#
#
# # Declare data and files location
# conf_file_address = '../muse_greenpeas.ini'
# obsData = sr.loadConfData(conf_file_address, group_variables=False)
# objList = obsData['sample_data']['object_list']
# fileList = obsData['sample_data']['file_list']
# dataFolder = obsData['sample_data']['data_folder']
# z_objs = obsData['sample_data']['z_array']
# voxel_list = [(170, 170), (92, 102)]
# mask_columns = ['wavelength', 'ion', 'w1', 'w2', 'w3', 'w4', 'w5', 'w6']
#
# for i, obj in enumerate(objList):
#
#     # Load the data
#     print(f'\n- {obj}')
#     fits_address_i = f'{dataFolder}/{fileList[i]}'
#     mask_address_i = f'{dataFolder}/{fileList[i].replace(".fits", "_mask.txt")}'
#     wave, cube, header = sr.import_fits_data(fits_address_i, instrument='MUSE')
#     wave_rest = wave / (1 + z_objs[i])
#
#     # Declare voxels
#     idx_voxel = voxel_list[i]
#     voxel = cube[:, idx_voxel[0], idx_voxel[1]]
#     flux_voxel = voxel.data.data
#
#     lm = sr.LineMesurer(wave_rest, flux_voxel)
#     lm.plot_spectrum_components()

    # # Identify the emission lines
    # noise_region = obsData['sample_data']['noiseRegion_array']
    # norm_flux = lm.continuum_remover(noise_region)
    # obsLinesTable = lm.line_finder(norm_flux, noiseWaveLim=noise_region, intLineThreshold=3)
    # obsLinesDF = lm.match_lines(obsLinesTable, sr._linesDb)
    # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=obsLinesDF)
    #
    # # Save the mask
    # idcsObsLines = (obsLinesDF.observation == 'detected')
    # lm.save_lineslog(obsLinesDF.loc[idcsObsLines, mask_columns], mask_address_i)
