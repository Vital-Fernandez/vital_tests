import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams
from astropy.wcs import WCS
from src.specsiser.print.plot import STANDARD_PLOT
from astro.data.muse.common_methods import compute_line_flux_image, image_array_binning

# Declare data and files location
obsData = sr.loadConfData('muse_J0925.ini', group_variables=False)
objList = np.array([obsData['sample_data']['object_list']])
fileList = np.array([obsData['sample_data']['file_list']])
dataFolder = Path(obsData['sample_data']['data_folder'])
resultsFolder = Path(obsData['sample_data']['results_folder'])
z_objs = np.array([obsData['sample_data']['z_array']])
pertil_array = obsData['sample_data']['percentil_array']
db_headers = obsData['sample_data']['database_header_list']
db_format = {'DEC_deg': '{: 0.8f}', 'RA_deg': '{: 0.8f}'}

lineAreas = {'H1_6563A_b': (6533.0, 6596.0),
             'S3_6312A': (6310.0, 6319.0),
             'O3_5007A': (4999.0, 5025.0),
             'O3_4363A': (4355.0, 4374.0)}

export_elements = ['intg_flux', 'intg_err', 'gauss_flux', 'gauss_err', 'amp', 'mu', 'sigma', 'amp_err', 'mu_err',
                   'sigma_err', 'v_r', 'v_r_err', 'sigma_vel', 'sigma_err_vel']

for i, obj in enumerate(objList):

    # Data location
    cube_address_i = dataFolder/fileList[i]
    objFolder = resultsFolder
    db_addresss = resultsFolder/f'{obj}_database.txt'

    # Load data
    obj_db = pd.read_csv(db_addresss, delim_whitespace=True, header=0, index_col=0)
    # wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')
    # wave = wave / (1 + z_objs[i])
    # print(f'\n- {obj}: Cube dimensions {cube.shape}')

    # Plot the line flux maps
    cube_shape = obsData['sample_data']['cube_size_array']
    for lineComp in obsData['default_line_fitting']['H1_6563A_b'].split('-'):
        for param in ['gauss_flux', 'sigma']:

            column_name = f'{lineComp}-{param}'
            lineFlux_i = np.reshape(obj_db[column_name].values, cube_shape.astype(int))

            # Define image countours based on the flux percentiles
            levelFlux_i = np.percentile(lineFlux_i[lineFlux_i > 0], pertil_array)
            levels_text_i = ['None'] * len(levelFlux_i)
            for idx, per in enumerate(pertil_array):
                levels_text_i[idx] = f'{levelFlux_i[idx]:.2f} $P_{{{per}}}$'

            # Plot line image map with coordinates
            labelsDict = {'xlabel': r'RA',
                          'ylabel': r'DEC',
                          'title': r'Galaxy {} {}'.format(obj, column_name)}

            # Plot Configuration
            defaultConf = STANDARD_PLOT.copy()
            defaultConf.update(labelsDict)
            rcParams.update({})

            # # Selecting plotting value pixels
            # frame_size = lineFlux_i.shape
            # x, y = np.arange(0, frame_size[1]), np.arange(0, frame_size[0])
            # X, Y = np.meshgrid(x, y)

            fig = plt.figure(figsize=(12, 8))
            idcs_non_nan = ~np.isnan(lineFlux_i)
            # ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
            ax = fig.add_subplot()
            im = ax.imshow(lineFlux_i)
            ax.set_title(column_name)
            fig.colorbar(im)
            plt.show()

            # fig = plt.figure(figsize=(12, 8))
            # # ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
            # ax = fig.add_subplot()
            #
            # CS3 = ax.contourf(X, Y, lineFlux_i, levels=levelFlux_i)
            # cbar = fig.colorbar(CS3)
            # cbar.ax.set_yticklabels(levels_text_i)
            # ax.set_facecolor('black')
            # ax.update(labelsDict)
            # imageName = f'{obj}_{column_name}_contours.png'
            # # plt.savefig(objFolder/imageName, bbox_inches='tight')
            # plt.show()

