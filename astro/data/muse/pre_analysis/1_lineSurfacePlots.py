import numpy as np
import src.specsiser as sr
from pathlib import Path
from mpdaf.obj import Cube
import matplotlib.pyplot as plt
import astropy.units as u
from mpdaf.obj import deg2sexa
from astropy.wcs import WCS


# Declare data and files location
conf_file_address = '../muse_greenpeas.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)
objList = obsData['sample_data']['object_list']
fileList = obsData['sample_data']['file_list']
dataFolder = obsData['sample_data']['data_folder']
idx_voxel = (170, 170)
lineAreas = {'H1_6563A': (6558.0, 6568.0), 'O3_5007A': (5002.0, 5013.0), 'S2_6717A': (6717.0, 6734.0)}
z_objs = obsData['sample_data']['z_array']


for i, obj in enumerate(objList):

    # Load the data
    print(f'\n- {obj}')
    file_address_i = f'{dataFolder}/{fileList[i]}'
    wave, cube, header = sr.import_fits_data(file_address_i, instrument='MUSE')
    wave = wave / (1 + z_objs[i])
    print(f'\n-- {header["OBJECT"]}')

    # Get astronomical coordinates one pixel
    coord_sky = cube.wcs.pix2sky(idx_voxel, unit=u.deg)
    dec, ra = deg2sexa(coord_sky)[0]
    wcs_cube = WCS(cube.data_header)

    # Treat all the lines
    for lineLabel, lineLimits in lineAreas.items():

        # Get line flux region
        ion, lineWave, latexlabel = sr.label_decomposition([lineLabel])
        lineIdcs = np.searchsorted(wave, np.array(lineLimits))
        lineImage = cube[lineIdcs[0]:lineIdcs[1], :, :].sum(axis=0)
        flux_image = lineImage.data.data

        # Plot line image map with coordinates
        labelsDict = {'xlabel': r'RA',
                      'ylabel': r'DEC',
                      'title': r'Galaxy {} {}'.format(obj, latexlabel[0])}
        outputFile = file_address_i.replace('.fits', f'_{lineLabel}_map.png')
        sr.plot.image_frame(flux_image, wcs=wcs_cube, axes_conf=labelsDict, output_file=outputFile)

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
