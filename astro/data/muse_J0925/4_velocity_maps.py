import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams
from astropy.wcs import WCS
from astro.data.muse.common_methods import compute_line_flux_image, image_array_binning

STANDARD_PLOT = {'figure.figsize': (20, 14), 'axes.titlesize': 20, 'axes.labelsize': 20, 'legend.fontsize': 12,
                 'xtick.labelsize': 12, 'ytick.labelsize': 12}

param_label = dict(gauss_flux=r'$F(\lambda)_{gauss}\,(erg\,cm^{-2} s^{-1} \AA^{-1})$',
                   mu='$\mu (\AA)$',
                   v_r='$v_{r}(km/s)$',
                   sigma='$\sigma(\AA)$',
                   sigma_vel='$\sigma(km/s)$')


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

for i, obj in enumerate(objList):

    # Data location
    cube_address_i = dataFolder/fileList[i]
    objFolder = resultsFolder
    db_addresss = resultsFolder/f'{obj}_database.txt'

    # Load data
    obj_db = pd.read_csv(db_addresss, delim_whitespace=True, header=0, index_col=0)
    wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')

    # Plot the line flux maps
    cube_shape = obsData['sample_data']['cube_size_array']
    for lineComp in obsData['default_line_fitting']['H1_6563A_b'].split('-'):
        for param in ['gauss_flux', 'mu', 'v_r', 'sigma', 'sigma_vel']:

            column_name = f'{lineComp}-{param}'
            lineFlux_i = np.reshape(obj_db[column_name].values, cube_shape.astype(int))
            ion, wave, latexCode = sr.label_decomposition(lineComp, scalar_output=True)

            # Define image countours based on the flux percentiles
            levelFlux_i = np.percentile(lineFlux_i[lineFlux_i > 0], pertil_array)
            levels_text_i = ['None'] * len(levelFlux_i)
            for idx, per in enumerate(pertil_array):
                levels_text_i[idx] = f'{levelFlux_i[idx]:.2f} $P_{{{per}}}$'

            # Crop image to the biggest square with non-nans
            nans = np.isnan(lineFlux_i)
            nancols = np.all(nans, axis=0)
            nanrows = np.all(nans, axis=1)
            firstcol = nancols.argmin()
            firstrow = nanrows.argmin()
            lastcol = len(nancols) - nancols[::-1].argmin()
            lastrow = len(nanrows) - nanrows[::-1].argmin()
            lineFlux_i_crop = lineFlux_i[firstrow:lastrow, firstcol:lastcol]

            # Image labels
            labelsDict = {'xlabel': r'RA',
                          'ylabel': r'DEC',
                          'title': r'Galaxy {}: {}'.format(obj, latexCode)}

            # Plot image
            rcParams.update(STANDARD_PLOT)

            fig = plt.figure(figsize=(12, 8))
            ax = fig.add_subplot()
            # ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))
            im = ax.imshow(lineFlux_i_crop)
            ax.set_title(column_name)
            ax.update(labelsDict)
            cbar = fig.colorbar(im)
            cbar.set_label(param_label[param], rotation=270,  size=16, labelpad=25)
            # cbar.ax.tick_params(labelsize=15)
            plt.show()
