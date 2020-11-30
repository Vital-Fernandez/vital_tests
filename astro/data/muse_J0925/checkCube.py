import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams
from astropy.wcs import WCS
from src.specsiser.print.plot import STANDARD_PLOT
from astro.data.muse.common_methods import compute_line_flux_image, image_array_binning, VoxelPlotter

lineAreas = {'H1_6563A': (6533.0, 6596.0),
             'S3_6312A': (6310.0, 6319.0),
             'O3_5007A': (4999.0, 5025.0),
             'O3_4363A': (4355.0, 4374.0)}

# Declare data and files location
obsData = sr.loadConfData('muse_J0925.ini', group_variables=False)
objList = np.array([obsData['sample_data']['object_list']])
fileList = np.array([obsData['sample_data']['file_list']])
dataFolder = Path(obsData['sample_data']['data_folder'])
resultsFolder = Path(obsData['sample_data']['results_folder'])
z_objs = np.array([obsData['sample_data']['z_array']])
pertil_array = obsData['sample_data']['percentil_array']
norm_flux = 1

for i, obj in enumerate(objList):

    # Data location
    cube_address_i = dataFolder/fileList[i]
    objFolder = resultsFolder/obj
    db_addresss = objFolder/f'{obj}_database.txt'

    # Load data
    wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')
    print(f'\n- {obj}: Cube dimensions {cube.shape}')

    # Get line region data
    lineFlux_dict, levelFlux_dict, levelText_dict = compute_line_flux_image(lineAreas,
                                                                            cube,
                                                                            z_objs[i],
                                                                            percent_array=pertil_array)

    # # -------- Plot the line flux maps
    # for lineLabel, lineLimits in lineAreas.items():
    #
    #     lineFlux_i = lineFlux_dict[lineLabel]
    #     levelFlux_i = levelFlux_dict[lineLabel]
    #     levelText_i = levelText_dict[lineLabel]
    #     flux_bin_i = image_array_binning(lineFlux_i, pertil_array)
    #
    #     # Plot line image map with coordinates
    #     labelsDict = {'xlabel': r'RA',
    #                   'ylabel': r'DEC',
    #                   'title': r'Galaxy {} {}'.format(obj, lineLabel)}
    #
    #     # Plot Configuration
    #     defaultConf = STANDARD_PLOT.copy()
    #     defaultConf.update(labelsDict)
    #     rcParams.update({})
    #
    #     # Selecting plotting value pixels
    #     frame_size = lineFlux_i.shape
    #     x, y = np.arange(0, frame_size[1]), np.arange(0, frame_size[0])
    #     X, Y = np.meshgrid(x, y)
    #
    #     fig = plt.figure(figsize=(12, 8))
    #     #ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
    #     ax = fig.add_subplot()
    #
    #     CS3 = ax.contourf(X, Y, lineFlux_i, levels=levelFlux_i)
    #     cbar = fig.colorbar(CS3)
    #     cbar.ax.set_yticklabels(levelText_i)
    #     ax.set_facecolor('black')
    #     ax.update(labelsDict)
    #     imageName = f'{obj}_{lineLabel}_contours.png'
    #     plt.show()
    #     # plt.savefig(objFolder/imageName, bbox_inches='tight')

    # -------- Plot a voxel
    voxel_coord = (149, 155)
    print(f'-- Treating voxel {voxel_coord}')
    idx_j, idx_i = voxel_coord

    # idx_database = (obj_db.y_voxel == idx_j) & (obj_db.x_voxel == idx_i)
    flux_voxel = cube[:, idx_j, idx_i].data.data * norm_flux
    flux_err = cube[:, idx_j, idx_i].var.data * norm_flux

    lm = sr.LineMesurer(wave, flux_voxel, redshift=z_objs[i], normFlux=norm_flux)
    lm.plot_spectrum_components()

    # -------- Plot interactive spectrum
    plotConf = {'image': {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'Galaxy {obj} O3_5007A'}}

    plotter = VoxelPlotter(wave, cube, lineFlux_dict['H1_6563A'], voxel_coord, image_fg=lineFlux_dict['O3_5007A'],
                           flux_levels=levelFlux_dict['O3_5007A'][3:], ax_user_conf=plotConf)