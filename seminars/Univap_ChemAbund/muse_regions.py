import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from src.specsiser.print.plot import STANDARD_PLOT
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, background_color, DARK_PLOT
from astropy.io import fits

folder = Path('/home/vital/Dropbox/Astrophysics/Seminars/UniVapo 2021/')

# Declare data and files location
conf_file = Path('/home/vital/PycharmProjects/vital_tests/astro/data/muse/muse_greenpeas.ini')
obsData = sr.loadConfData(conf_file, group_variables=False)
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
defaultConf = DARK_PLOT.copy()
rcParams.update(defaultConf)

for i, obj in enumerate(objList):

    # Data location
    cube_address_i = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    db_addresss = objFolder/f'{obj}_database.fits'

    # Load data
    wave, cube, header = sr.import_fits_data(cube_address_i, instrument='MUSE')
    wave = wave / (1 + z_objs[i])
    print(f'\n- {obj}: Cube dimensions {cube.shape}')

    # Create flux maps for the main lines:
    for lineLabel, lineLimits in lineAreas.items():

        plot_image_file = objFolder/f'{obj}_{lineLabel}_contours.png'
        ion, wavelength, latexLabel = sr.label_decomposition(lineLabel, scalar_output=True)

        # Extract cube slice using mpdaf defult tools.
        # This requires the input wavelengths to be on the same scale as in the cube
        line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
        flux_image = line_image.data.data
        levelContours = np.nanpercentile(flux_image, pertil_array)

        # loop throught the intensity layers
        # for idx_ion_boundary, flux_perdentil(levelContours):
        idx_ion_boundary=2

        # Plot the image:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))

        im = ax.imshow(flux_image, cmap=cm.gray, norm=colors.SymLogNorm(linthresh=levelContours[-3], vmin=levelContours[-3], base=10))
        cntr1 = ax.contour(flux_image, levels=levelContours[::-1][idx_ion_boundary:], cmap='viridis', norm=colors.LogNorm())

        for idx, percentile in enumerate(pertil_array[::-1][idx_ion_boundary:]):
            label = r'$P_{{{}}}([SIII]6312\AA)$'.format(idx, percentile)
            cntr1.collections[idx].set_label(label)
        ax.legend()

        ax.update({'title': r'{} galaxy, {} flux'.format(obj, latexLabel), 'xlabel': r'RA', 'ylabel': r'DEC'})
        # plt.savefig(plot_image_file, bbox_inches='tight')
        plt.tight_layout
        plt.show()

