import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams
from astropy.wcs import WCS
from src.specsiser.print.plot import STANDARD_PLOT
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, image_array_binning
from astropy.io import fits

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
coordinates_keys_list = obsData['data_location']['wcs_key_list']

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']

# Plot set up
labelsDict = {'xlabel': r'RA',
              'ylabel': r'DEC'}
defaultConf = STANDARD_PLOT.copy()
defaultConf.update(labelsDict)
rcParams.update({})

for i, obj in enumerate(objList):

    # Data location
    cube_address_i = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    db_addresss = objFolder/f'{obj}_database.fits'

    # Load data
    wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')
    wave = wave / (1 + z_objs[i])
    print(f'\n- {obj}: Cube dimensions {cube.shape}')

    # Create empty fits file
    new_hdul = fits.HDUList()
    new_hdul.append(fits.PrimaryHDU())

    # Second page for the fits file plot configuration
    col_waves = fits.Column(name='wave', array=wave, format='1E')
    hdu_table = fits.BinTableHDU.from_columns([col_waves], name='PlotConf')
    new_hdul.append(hdu_table)
    for key in coordinates_keys_list:
        new_hdul[1].header[key] = cube.data_header[key]
    new_hdul[1].header['NPIXWAVE'] = cube.data_header['NAXIS3']

    # Create flux maps for the main lines:
    for lineLabel, lineLimits in lineAreas.items():

        # Extract cube slice using mpdaf defult tools.
        # This requires the input wavelengths to be on the same scale as in the cube
        line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
        flux_image = line_image.data.data

        # Personal scale for the image flux
        log_flux = np.log10(flux_image)
        flux_contours = np.zeros(flux_image.shape)
        idcs_removed = np.logical_or(log_flux < 0.0, np.isnan(log_flux))
        flux_contours[~idcs_removed] = log_flux[~idcs_removed]

        # Store fluxes and contours
        new_hdul.append(fits.ImageHDU(name=f'{lineLabel}_flux', data=flux_image, ver=1))
        new_hdul.append(fits.ImageHDU(name=f'{lineLabel}_contour', data=flux_contours, ver=1))

        # Define image countours based on the flux percentiles
        levelFlux_i = np.percentile(flux_contours[flux_contours > 0], pertil_array)
        levels_text_i = ['None'] * len(levelFlux_i)
        for idx, per in enumerate(pertil_array):
            levels_text_i[idx] = f'{levelFlux_i[idx]:.2f} $P_{{{per}}}$ $log(F_{{\lambda}})$'

        # Plot the image:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))

        frame_size = flux_image.shape
        x, y = np.arange(0, frame_size[1]), np.arange(0, frame_size[0])
        X, Y = np.meshgrid(x, y)

        CS3 = ax.contourf(X, Y, flux_contours, levels=levelFlux_i)
        cbar = fig.colorbar(CS3)
        cbar.ax.set_yticklabels(levels_text_i)
        ax.set_facecolor('black')
        ax.update(labelsDict)
        ax.set_title(f'Galaxy {obj} {lineLabel}')
        imageName = f'{obj}_{lineLabel}_contours.png'
        plt.savefig(objFolder/imageName, bbox_inches='tight')

    # Store the drive
    new_hdul.writeto(db_addresss, overwrite=True, output_verify='fix')

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
# idx_voxel = (170, 170)
# lineAreas = {'H1_6563A': (6558.0, 6568.0), 'O3_5007A': (5002.0, 5013.0), 'S2_6717A': (6717.0, 6734.0)}
# z_objs = obsData['sample_data']['z_array']
#
#
# for i, obj in enumerate(objList):
#
#     # Load the data
#     print(f'\n- {obj}')
#     file_address_i = f'{dataFolder}/{fileList[i]}'
#     wave, cube, header = sr.import_fits_data(file_address_i, instrument='MUSE')
#     wave = wave / (1 + z_objs[i])
#     print(f'\n-- {header["OBJECT"]}')
#
#     # Get astronomical coordinates one pixel
#     coord_sky = cube.wcs.pix2sky(idx_voxel, unit=u.deg)
#     dec, ra = deg2sexa(coord_sky)[0]
#     wcs_cube = WCS(cube.data_header)
#
#     # Treat all the lines
#     for lineLabel, lineLimits in lineAreas.items():
#
#         # Get line flux region
#         ion, lineWave, latexlabel = sr.label_decomposition([lineLabel])
#         lineIdcs = np.searchsorted(wave, np.array(lineLimits))
#         lineImage = cube[lineIdcs[0]:lineIdcs[1], :, :].sum(axis=0)
#         flux_image = lineImage.data.data
#
#         # Plot line image map with coordinates
#         labelsDict = {'xlabel': r'RA',
#                       'ylabel': r'DEC',
#                       'title': r'Galaxy {} {}'.format(obj, latexlabel[0])}
#         outputFile = file_address_i.replace('.fits', f'_{lineLabel}_map.png')
#         sr.plot.image_frame(flux_image, wcs=wcs_cube, axes_conf=labelsDict)#, output_file=outputFile)

# # Load data
# files_list = ['CGCG007.fits', 'UGC5205.fits']
# files_path = Path('D:/Google drive/Astrophysics/Datos/MUSE - Amorin')
# files_address = list(files_path/file for file in files_list)
# cube = Cube(filename=str(files_address[0]))
#
# # Cube shape (lambda, Y, X) // Header shape (1, 2, 3) = (RA---TAN, DEC--TAN, AWAV) = (X, Y, lambda) // QfitsView (X, Y)
# cube_size = cube.shape
#
# # Fits properties
# cube.info()
#
# # Fits header
# hdr = cube.data_header
#
# # Reconstruct wave:
# cube.wave.info()
# dw = hdr['CD3_3']
# w_min = hdr['CRVAL3']
# nPixels = hdr['NAXIS3']
# w_max = w_min + dw * nPixels
# wave = np.linspace(w_min, w_max, nPixels, endpoint=False)
#
# # Redshift correction
# z = 0.0046
# wave = wave / (1 + z)
#
# # Get voxel spectrum
# idx_voxel = (170, 170)
# voxel = cube[:, idx_voxel[0], idx_voxel[1]]
# flux_voxel = voxel.data.data
#
# # Get line flux region
# lineLabel, lineWave = 'H1_6563A', 6563.0
# lineRegions = np.array([lineWave-5, lineWave, lineWave+5])
# lineIdcs = np.searchsorted(wave, lineRegions)
# lineImage = cube[lineIdcs[0]:lineIdcs[2], :, :].sum(axis=0)
# flux_image = lineImage.data.data
#
# # Get astronomical coordinates one pixel
# coord_sky = cube.wcs.pix2sky(idx_voxel, unit=u.deg)
# dec, ra = deg2sexa(coord_sky)[0]
# wcs_cube = WCS(cube.data_header)
#
# # Plot the data
# labelsDict = {'xlabel': r'Wavelength $(\AA)$',
#               'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1}) \cdot 10^{20}$',
#               'title': f'Galaxy CGCG007 (voxel coords: {idx_voxel[0]}, {idx_voxel[1]})'}
# sr.plot.spectrum(wave, flux_voxel, axes_conf=labelsDict)
#
# # Plot line image map
# labelsDict = {'xlabel': r'X',
#               'ylabel': r'Y',
#               'title': r'Galaxy CGCG007 $H\alpha$'}
# sr.plot.image_frame(flux_image, axes_conf=labelsDict)
#
# # Plot line image map with coordinates
# labelsDict = {'xlabel': r'RA',
#               'ylabel': r'DEC',
#               'title': r'Galaxy CGCG007 $H\alpha$'}
# sr.plot.image_frame(flux_image, wcs=wcs_cube, axes_conf=labelsDict)
#
# # Plot line image contours
# labelsDict = {'xlabel': r'RA',
#               'ylabel': r'DEC',
#               'title': r'Galaxy CGCG007 $H\alpha$'}
# sr.plot.image_contour(flux_image, wcs=wcs_cube, axes_conf=labelsDict)




# coord_sky = cube.wcs.pix2sky([8, 28], unit=u.deg)
# dec, ra = deg2sexa(coord_sky)[0]
# voxel = cube[:, idx_voxel[0], idx_voxel[1]]
# voxel.plot(title = 'Zoom Spectrum ra=%s dec=%s' %(ra, dec))
# plt.show()
# sr.plot.spectrum()
# spectrum.fit_lines()

# print(cube.shape)
# print(cube.data.shape)
# print(cube.var.shape)
# print(cube.mask.shape)
