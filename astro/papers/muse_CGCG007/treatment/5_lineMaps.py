import lime
from pathlib import Path
from astropy.io import fits
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits, param_images

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
coordinates_keys_list = obsData['data_location']['wcs_key_list']

# Masks to consider
masks = 'all'

for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    voxelFolder = resultsFolder/obj/'voxel_data'
    db_address = objFolder / f'{obj}_database.fits'
    fitsLog_address = objFolder / f'{obj}_linesLog.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Load data
    wave, cube, header = import_muse_fits(cube_address)

    # Header data with the cosmological coordinates
    hdr_plot = fits.getheader(db_address, extname='PlotConf')

    plot_dict = {'CRPIX1': header['CRPIX1'],
                 'CRPIX2': header['CRPIX2'],
                 'CD1_1': header['CD1_1'],
                 'CD1_2': header['CD1_2'],
                 'CD2_1': header['CD2_1'],
                 'CD2_2': header['CD2_2'],
                 'CUNIT1': header['CUNIT1'],
                 'CUNIT2': header['CUNIT2'],
                 'CTYPE1': header['CTYPE1'],
                 'CTYPE2': header['CTYPE2']}

    # Generate the maps
    lime.save_param_maps(fitsLog_address, param_images, objFolder, maskFits_address, ext_mask=masks, ext_log='_LINELOG',
                         page_hdr=plot_dict)


    # param = 'gauss_flux'
    # user_lines = param_images[param]
    #
    # fits_file = Path(objFolder) / f'{param}.fits'
    # with fits.open(fits_file):
    #     for line in user_lines:
    #         param_image = fits.getdata(fits_file, line)
    #         param_hdr = fits.getheader(fits_file, line)
    #
    #         fig = plt.figure(figsize=(10, 10))
    #         ax = fig.add_subplot(projection=WCS(fits.Header(param_hdr)), slices=('x', 'y'))
    #         im = ax.imshow(param_image)
    #         ax.update({'title': f'Galaxy {obj}: {param}-{line}', 'xlabel': r'RA', 'ylabel': r'DEC'})
    #         plt.show()