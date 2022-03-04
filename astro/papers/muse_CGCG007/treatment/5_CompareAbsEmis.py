import numpy as np
import lime
from pathlib import Path
from astropy.io import fits
from astro.papers.muse_CGCG007.muse_CGCG007_methods import import_fado_cube, import_muse_fits

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList_obs = obsData['data_location']['file_list']
fileList_abs = obsData['data_location']['file_abs_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
fitsFolder_abs = Path(obsData['data_location']['fits_abs_folder'])

z_objs = obsData['sample_data']['z_array']
norm_flux_abs = obsData['sample_data']['norm_flux_abs']
norm_flux = obsData['sample_data']['norm_flux']

for i, obj in enumerate(objList):

    # Data location
    cube_address_obs = fitsFolder/fileList_obs[i]
    cube_address_abs = fitsFolder_abs/fileList_abs[i]
    objFolder = resultsFolder/obj
    ObsLog_address = objFolder / f'{obj}_linesLog.fits'
    AbsLog_address = objFolder / f'{obj}_linesLog_abs.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    idxY, idxX = 167, 167
    obs_log = lime.load_lines_log(ObsLog_address, ext=f'{idxY}-{idxX}_LINELOG')
    abs_log = lime.load_lines_log(AbsLog_address, ext=f'{idxY}-{idxX}_LINELOG')

    for line in ['H1_6563A', 'H1_4861A', 'He1_5016A', 'H1_9015A']:

        norm = obs_log.loc[line, "cont"]/abs_log.loc[line, "cont"]

        print('\n', line, f'{norm}')
        print(f'- Gauss Flux : Emis =  {obs_log.loc[line, "gauss_flux"]}, Abs = {abs_log.loc[line, "gauss_flux"] * norm}')
        print(f'- Cont Flux :  Emis =  {obs_log.loc[line, "cont"]}, Abs = {abs_log.loc[line, "cont"] * norm}')

        1 - (obs_log.loc[line, "gauss_flux"]/(-abs_log.loc[line, "gauss_flux"]*norm)) * 100

    param = 'gauss_flux'
    user_lines = param_images[param]

    fits_file = Path(objFolder) / f'{param}.fits'
    with fits.open(fits_file):
        for line in user_lines:
            param_image = fits.getdata(fits_file, line)
            param_hdr = fits.getheader(fits_file, line)

            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(projection=WCS(fits.Header(param_hdr)), slices=('x', 'y'))
            im = ax.imshow(param_image)
            ax.update({'title': f'Galaxy {obj}: {param}-{line}', 'xlabel': r'RA', 'ylabel': r'DEC'})
            plt.show()