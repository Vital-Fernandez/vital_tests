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
from astropy.io import fits
from src.specsiser.print.plot import STANDARD_PLOT
from astropy.visualization import mpl_normalize, SqrtStretch
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, VoxelPlotter, reconstruct_wavelength


# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']
dict_errs = {}

fits_folder = Path('/home/vital/Astro-data/Observations/MUSE - Amorin/FADO_analysis/Z4SalpP2000/')
fits_3 = 'cgcg007025_HBIN024_FDres2_3DnebSED.fits'
fits_4 = 'cgcg007025_HBIN024_FDres2_3DnoNEB.fits'
fits_5 = 'cgcg007025_HBIN024_FDres2_3DOBS.fits'
fits_6 = 'cgcg007025_HBIN024_FDres2_3DstelFIT.fits'
file_address, ext = fits_folder/fits_5, 0

with fits.open(file_address) as hdu_list:
    data_fado = hdu_list[ext].data
    hdr = hdu_list[ext].header
wave_fado = reconstruct_wavelength(hdr)

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        db_address_i = dataFolder/f'{obj}_database.txt'
        cube_address_i = fitsFolder/fileList[i]
        output_folder = dataFolder/obj
        mask_global_address_i = dataFolder / f'{obj}_mask.txt'

        # Load data
        wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')
        wave_rest = wave / (1 + z_objs[i])

        # Reset and measure the lines
        # lm = sr.LineMesurer(wave_rest, flux_voxel * norm_flux, normFlux=norm_flux, linesDF_address=mask_global_address_i)

        print(f'\n- {obj}: Cube dimensions {cube.shape}')

        # Get line region data
        lineFlux_dict, levelFlux_dict, levelText_dict = compute_line_flux_image(lineAreas, cube, z_objs[i],
                                                                                percent_array=pertil_array)
        # Declare voxels to analyse
        lineLabel = 'O3_5007A'
        percentil_array = pertil_array
        wcs_cube = WCS(cube.data_header)
        voxel_coord = (173, 169)
        idx_j, idx_i = voxel_coord

        fluxImage = lineFlux_dict[lineLabel]
        fluxLevels = levelFlux_dict[lineLabel]
        levels_text = levelText_dict[lineLabel]
        flux_voxel = cube[:, idx_j, idx_i].data.data

        plotConf = {'image': {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'Galaxy {obj} {lineLabel}'}}

        plotter = VoxelPlotter(wave_fado, data_fado, lineFlux_dict['H1_6563A'], voxel_coord, image_fg=lineFlux_dict[lineLabel],
                              flux_levels=fluxLevels[::-1][2:], ax_user_conf=plotConf, header=cube.data_header)
        # plotter.plot_map_voxel(lineFlux_dict['H1_6563A'], voxel_coord, image_fg=lineFlux_dict[lineLabel],
        #                        flux_levels=fluxLevels[2:])
        # plt.show()

