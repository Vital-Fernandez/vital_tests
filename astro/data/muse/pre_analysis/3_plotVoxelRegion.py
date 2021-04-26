import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import src.specsiser as sr
from pathlib import Path
from mpdaf.obj import Cube
from matplotlib import pyplot as plt, rcParams
import astropy.units as u
from mpdaf.obj import deg2sexa
from astropy.wcs import WCS
from src.specsiser.print.plot import STANDARD_PLOT
from astropy.visualization import mpl_normalize, SqrtStretch
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, VoxelPlotter

from matplotlib import pyplot as plt, rcParams, gridspec
from src.specsiser.physical_model.line_tools import STANDARD_PLOT, STANDARD_AXES



# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['sample_data']['object_list']
fileList = obsData['sample_data']['file_list']
dataFolder = Path(obsData['sample_data']['data_folder'])
z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']
dict_errs = {}

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        db_address_i = dataFolder/f'{obj}_database.txt'
        cube_address_i = dataFolder/fileList[i]
        output_folder = dataFolder/obj
        mask_global_address_i = dataFolder / f'{obj}_mask.txt'

        # Load data
        wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')
        wave_rest = wave / (1 + z_objs[i])
        obj_db = pd.read_csv(db_address_i, delim_whitespace=True, header=0, index_col=0)
        mask_df = pd.read_csv(mask_global_address_i, delim_whitespace=True, header=0, index_col=0)

        # Reset and measure the lines
        # lm = sr.LineMesurer(wave_rest, flux_voxel * norm_flux, normFlux=norm_flux, linesDF_address=mask_global_address_i)

        print(f'\n- {obj}: Cube dimensions {cube.shape}')

        # Get line region data
        lineFlux_dict, levelFlux_dict, levelText_dict = compute_line_flux_image(lineAreas, cube, z_objs[i],
                                                                                percent_array=pertil_array)
        # Declare voxels to analyse
        lineLabel = 'S3_6312A'
        percentil_array = pertil_array
        wcs_cube = WCS(cube.data_header)
        voxel_coord = (173, 169)
        idx_j, idx_i = voxel_coord

        fluxImage = lineFlux_dict[lineLabel]
        fluxLevels = levelFlux_dict[lineLabel]
        levels_text = levelText_dict[lineLabel]
        flux_voxel = cube[:, idx_j, idx_i].data.data

        plotConf = {'image': {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'Galaxy {obj} {lineLabel}'}}

        plotter = VoxelPlotter(wave_rest, cube, lineFlux_dict['H1_6563A'], voxel_coord, image_fg=lineFlux_dict[lineLabel],
                              flux_levels=fluxLevels[2:], ax_user_conf=plotConf)
        # plotter.plot_map_voxel(lineFlux_dict['H1_6563A'], voxel_coord, image_fg=lineFlux_dict[lineLabel],
        #                        flux_levels=fluxLevels[2:])
        # plt.show()

