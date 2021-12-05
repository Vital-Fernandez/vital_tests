import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astropy.io import fits
import lime as lm
from astro.papers.SHOC579_project.SHOC579_methods import open_manga_cubes, line_regions, STANDARD_PLOT

# # Declare data and files location
# obsData = sr.loadConfData('../muse_CGCG007.ini')
# objList = obsData['data_location']['object_list']
# fileList = obsData['data_location']['file_list']
# fitsFolder = Path(obsData['data_location']['fits_folder'])
# dataFolder = Path(obsData['data_location']['data_folder'])
# resultsFolder = Path(obsData['data_location']['results_folder'])

obs_conf = lm.load_cfg(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\obsConf.ini')

files_data = obs_conf['data_location_windows']
data_folder = Path(files_data['data_folder'])
results_folder = Path(files_data['results_folder'])
file_list = files_data['file_list']
obj_list = files_data['obj_list']

z_objs = obs_conf['sample_data']['z_array']
percentil_array = obs_conf['sample_data']['percentil_array']
coordinates_keys_list = obs_conf['sample_data']['wcs_key_list']

# Plot set up
labelsDict = {'xlabel': r'RA',
              'ylabel': r'DEC'}
defaultConf = STANDARD_PLOT.copy()
defaultConf.update(labelsDict)
rcParams.update({})

for i, obj in enumerate(obj_list):

    # Data location
    cube_address_i = data_folder/file_list[i]
    objFolder = results_folder/obj
    db_address = objFolder/f'{obj}_database.fits'

    # Load data
    wave, flux, err, hdr = open_manga_cubes(cube_address_i)
    print(f'\n- {obj}: Cube dimensions {flux.shape}')

    # Create empty fits file
    new_hdul = fits.HDUList()
    new_hdul.append(fits.PrimaryHDU())

    # Second page for the fits file plot configuration
    col_waves = fits.Column(name='wave', array=wave, format='1E')
    hdu_table = fits.BinTableHDU.from_columns([col_waves], name='PlotConf')
    new_hdul.append(hdu_table)
    for key in coordinates_keys_list:
        if key in hdr:
            new_hdul[1].header[key] = hdr[key]
    new_hdul[1].header['NPIXWAVE'] = hdr['NAXIS3']

    # Get line regions data
    idcs_Halpha = np.searchsorted(wave, line_regions['H1_6563A'])
    Halpha_flux = flux[idcs_Halpha[0]:idcs_Halpha[1], :, :].sum(axis=0)
    levelContoursHalpha = np.nanpercentile(Halpha_flux, percentil_array)

    # Create flux maps for the main lines:
    for lineLabel, lineLimits in line_regions.items():

        plot_image_file = objFolder/f'{obj}_{lineLabel}_contours.png'
        ion, wavelength, latexLabel = lm.label_decomposition(lineLabel, scalar_output=True)

        # Extract cube slice using mpdaf defult tools.
        idcs_line = np.searchsorted(wave, lineLimits)
        lineMap = flux[idcs_line[0]:idcs_line[1], :, :].sum(axis=0)
        levelContours = np.nanpercentile(lineMap, percentil_array)

        # Store fluxes and contours
        hdu_image = fits.ImageHDU(name=f'{lineLabel}_flux', data=lineMap, ver=1)
        for idx, level in enumerate(levelContours):
            level_label = f'hierarch P{int(percentil_array[idx]*100)}'
            hdu_image.header[level_label] = level
        new_hdul.append(hdu_image)

        # Plot the image:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y', 1))

        min_background_percentil = levelContoursHalpha[1]
        normalization_background = colors.SymLogNorm(linthresh=min_background_percentil,
                                                     vmin=min_background_percentil,
                                                     base=10)
        im = ax.imshow(Halpha_flux, cmap=cm.gray, norm=normalization_background)

        # Plot contours
        cntr1 = ax.contour(lineMap, levels=levelContours, cmap='viridis', norm=colors.LogNorm())
        for idx, percentile in enumerate(percentil_array[1:]):
            label = r'$P_{{{}}}$({})'.format(percentil_array[idx], latexLabel)
            cntr1.collections[idx].set_label(label)
        ax.legend()

        ax.update({'title': r'{} galaxy, {} flux'.format(obj, latexLabel), 'xlabel': r'RA', 'ylabel': r'DEC'})
        plt.savefig(plot_image_file, bbox_inches='tight')
        # plt.show()

    # Store the drive
    new_hdul.writeto(db_address, overwrite=True, output_verify='fix')


