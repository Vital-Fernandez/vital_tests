import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors
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

    # Load data
    wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')
    wave = wave / (1 + z_objs[i])
    print(f'\n- {obj}: Cube dimensions {cube.shape}')

    lineLabel = 'H1_6563A'
    lineLimits = lineAreas[lineLabel]

    # This requires the input wavelengths to be on the same scale as in the cube
    line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
    flux_image = line_image.data.data

    levelContours = np.nanpercentile(flux_image, pertil_array)
    for i, percentile in enumerate(levelContours):
        print(f'{i}, Level P({pertil_array[i]}) = {percentile} flux')

    maFlux_image = np.ma.masked_where((flux_image >= levelContours[-5]) &
                                      (flux_image < levelContours[-4]),
                                      flux_image)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
    # im = ax.imshow(flux_image, vmin=0, interpolation='none')
    # im = ax.imshow(flux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=levelContours[1], vmin=levelContours[1], base=10))
    im = ax.imshow(maFlux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=levelContours[1], vmin=levelContours[1], base=10))
    # im = ax.imshow(flux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=levelContours[1], vmin=levelContours[2]))
    ax.contour(flux_image, levels=[levelContours[-5]], colors='red', alpha=0.5)

    # cbar = fig.colorbar(im, ticks=levelContours)
    print(levelContours[2:6])
    ax.set_title(f'Galaxy {obj} {lineLabel}')
    plt.show()


    # Image with a clear display of the galaxy features
    # lineLabel = 'H1_6563A'
    # lineLimits = lineAreas[lineLabel]
    #
    # # This requires the input wavelengths to be on the same scale as in the cube
    # line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
    # flux_image = line_image.data.data
    #
    #
    # levelContours = np.nanpercentile(flux_image, pertil_array)
    # for i, percentile in enumerate(levelContours):
    #     print(f'{i}, Level P({pertil_array[i]}) = {percentile} flux')
    #
    # fig = plt.figure(figsize=(12, 8))
    # ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
    # # im = ax.imshow(flux_image, vmin=0, interpolation='none')
    # im = ax.imshow(flux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=levelContours[1], vmin=levelContours[1], base=10))
    # # im = ax.imshow(flux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=levelContours[1], vmin=levelContours[2]))
    # ax.contour(flux_image, levels=[levelContours[-5]], colors='white', alpha=0.5)
    #
    # # cbar = fig.colorbar(im, ticks=levelContours)
    # print(levelContours[2:6])
    # ax.set_title(f'Galaxy {obj} {lineLabel}')
    # plt.show()


    # Image with the galaxy and the two contours exposition
    # lineLabel = 'H1_6563A'
    # lineLimits = lineAreas[lineLabel]
    #
    # # This requires the input wavelengths to be on the same scale as in the cube
    # line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
    # flux_image = line_image.data.data
    #
    # S3_6312A_image = cube.get_image(np.array(lineAreas['S3_6312A']) * (1 + z_objs[i]), subtract_off=True)
    # S3_6312A_flux = S3_6312A_image.data.data
    # S3_6312A_percentiles = np.nanpercentile(S3_6312A_flux, pertil_array)
    #
    # levelContours = np.nanpercentile(flux_image, pertil_array)
    # for i, percentile in enumerate(levelContours):
    #     print(f'{i}, Level P({pertil_array[i]}) = {percentile} flux')
    #
    # print(np.log10(levelContours))
    # print(np.nanpercentile(np.log10(flux_image), pertil_array))
    #
    # fig = plt.figure(figsize=(12, 8))
    # ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
    # im = ax.imshow(flux_image, vmin=0, interpolation='none')
    # ax.contour(flux_image, levels=[levelContours[-5]], colors='white', alpha=0.5)
    # ax.contour(S3_6312A_flux, levels=S3_6312A_percentiles[-4:-1], cmap='Pastel1')
    #
    # # cbar = fig.colorbar(im, ticks=levelContours)
    # print(levelContours[2:6])
    # ax.set_title(f'Galaxy {obj} {lineLabel}')
    # plt.show()