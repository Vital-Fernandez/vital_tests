import numpy as np
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, cm, colors
from astropy.wcs import WCS
from astropy.io import fits

# Declare data and files location
objList = ['CGCG007']
fileList = ['CGCG007.fits']
fitsFolder = Path('/home/vital/Astro-data/Observations/MUSE - Amorin')
dataFolder = Path('/home/vital/Dropbox/Astrophysics/Papers/muse_CGCG007/data')
resultsFolder = Path('/home/vital/Dropbox/Astrophysics/Papers/muse_CGCG007/treatment')
coordinates_keys_list = ['CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','CUNIT1','CUNIT2','CTYPE1','CTYPE2',
                         'CSYER1','CSYER2','CTYPE3','CUNIT3','CD3_3','CRPIX3','CRVAL3','CRDER3','CD1_3','CD2_3',
                         'CD3_1','CD3_2']

z_objs = np.array([0.004691, 0.005113])
pertil_array = np.array([0, 90.50, 92.50, 95.50, 97.50, 99.50, 99.90, 99.99])

# Plot set up
labelsDict = {'xlabel': r'RA',
              'ylabel': r'DEC'}

lineAreas = {'H1_6563A': (6558.0, 6568.0),
             'S3_6312A': (6308.15, 6317.25),
             'O3_5007A': (5002.0, 5013.0),
             'S2_6717A': (6717.0, 6734.0)}

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

    # First page has the default (primary) data
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

        plot_image_file = objFolder/f'{obj}_{lineLabel}_contours.png'
        ion, wavelength, latexLabel = sr.label_decomposition(lineLabel, scalar_output=True)

        # Extract cube slice using mpdaf defult pyPopstar.
        # This requires the input wavelengths to be on the same scale as in the cube
        line_image = cube.get_image(np.array(lineLimits) * (1 + z_objs[i]), subtract_off=True)
        flux_image = line_image.data.data
        levelContours = np.nanpercentile(flux_image, pertil_array)

        # Store fluxes and contours
        hdu_image = fits.ImageHDU(name=f'{lineLabel}_flux', data=flux_image, ver=1)
        for idx, level in enumerate(levelContours):
            level_label = f'hierarch P{int(pertil_array[idx]*100)}'
            hdu_image.header[level_label] = level
        new_hdul.append(hdu_image)

        # Plot the image:
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection=WCS(cube.data_header), slices=('x', 'y', 1))

        min_background_percentil = levelContours[2]
        normalization_background = colors.SymLogNorm(linthresh=min_background_percentil, vmin=min_background_percentil, base=10)
        im = ax.imshow(flux_image, cmap=cm.gray, norm=normalization_background)

        # Plot contours
        min_percentil_background = 1
        contours_levels = levelContours[min_percentil_background:]
        cntr1 = ax.contour(flux_image, levels=contours_levels, cmap='viridis', norm=colors.LogNorm())
        for idx, percentile in enumerate(pertil_array[min_percentil_background:]):
            label = r'$P_{{{}}}$({})'.format(pertil_array[idx], latexLabel)
            cntr1.collections[idx].set_label(label)
        ax.legend()

        ax.update({'title': r'{} galaxy, {} flux'.format(obj, latexLabel), 'xlabel': r'RA', 'ylabel': r'DEC'})
        # plt.savefig(plot_image_file, bbox_inches='tight')
        plt.show()

        # Store the drive
        new_hdul.writeto(db_addresss, overwrite=True, output_verify='fix')

    # Make plot from stored fits
    hdr = fits.getheader(db_addresss, extname='PlotConf')
    for lineLabel, lineLimits in lineAreas.items():
        flux_image = fits.getdata(db_addresss, f'{lineLabel}_flux', ver=1)
        flux_contours = fits.getdata(db_addresss, f'{lineLabel}_contour', ver=1)

        # Define image countours based on the flux percentiles
        levelFlux_i = np.percentile(flux_contours[flux_contours > 0], pertil_array)
        levels_text_i = ['None'] * len(levelFlux_i)
        for idx, per in enumerate(pertil_array):
            levels_text_i[idx] = f'{levelFlux_i[idx]:.2f} $P_{{{per}}}$ $log(F_{{\lambda}})$'

        # Plot the image:
        fig = plt.figure(figsize=(12, 8))
        ax = fig.add_subplot(projection=WCS(hdr), slices=('x', 'y', 1))

        frame_size = flux_image.shape
        x, y = np.arange(0, frame_size[1]), np.arange(0, frame_size[0])
        X, Y = np.meshgrid(x, y)

        CS3 = ax.contourf(X, Y, flux_contours, levels=levelFlux_i)
        cbar = fig.colorbar(CS3)
        cbar.ax.set_yticklabels(levels_text_i)
        ax.set_facecolor('black')
        ax.update(labelsDict)
        ax.set_title(f'Galaxy {obj} {lineLabel}')
        ax.set_xlabel('RA')
        ax.set_ylabel('DEC')
        plt.show()
