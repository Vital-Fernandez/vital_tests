import lime
from pathlib import Path
from astropy.io import fits
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_muse_fits, import_fado_cube, param_images_abs

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList_obs = obsData['data_location']['file_list']
fileList_abs = obsData['data_location']['file_abs_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
fitsFolder_abs = Path(obsData['data_location']['fits_abs_folder'])

# Masks to consider
masks = 'all'

for i, obj in enumerate(objList):

    # Data location
    cube_address_obs = fitsFolder/fileList_obs[i]
    cube_address_abs = fitsFolder_abs/fileList_abs[i]
    objFolder = resultsFolder/obj
    fitsLog_address = objFolder / f'{obj}_linesLog_abs.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    # Load data
    wave_abs, cube_abs, header_abs = import_fado_cube(cube_address_abs)

    # Open emission cube
    wave_obs, cube_obs, header_obs = import_muse_fits(cube_address_obs)

    plot_dict = {'CRPIX1': header_obs['CRPIX1'],
                 'CRPIX2': header_obs['CRPIX2'],
                 'CD1_1': header_obs['CD1_1'],
                 'CD1_2': header_obs['CD1_2'],
                 'CD2_1': header_obs['CD2_1'],
                 'CD2_2': header_obs['CD2_2'],
                 'CUNIT1': header_obs['CUNIT1'],
                 'CUNIT2': header_obs['CUNIT2'],
                 'CTYPE1': header_obs['CTYPE1'],
                 'CTYPE2': header_obs['CTYPE2']}

    # Generate the maps
    lime.save_param_maps(fitsLog_address, param_images_abs, objFolder, maskFits_address, ext_mask=masks,
                         ext_log='_LINELOG', page_hdr=plot_dict, output_files_prefix='abs_')


