import numpy as np

import lime
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS

from astro.papers.muse_CGCG007.muse_CGCG007_methods import label_Conver, abs_target_lines, param_images_abs, abs_target_lines

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
voxel_grid_size = obsData['sample_data']['grid_shape_array']
coordinates_keys_list = obsData['data_location']['wcs_key_list']
# ------ Emission/absroption
for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    db_address = objFolder / f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    emis_fits = objFolder / f'gauss_flux.fits'
    abs_fits = objFolder/'abs_gauss_flux.fits'
    Hbeta_emis = fits.getdata(emis_fits, 'H1_4861A')
    Hbeta_abs = fits.getdata(abs_fits, 'H1_4861A')

    hdr_plot = fits.getheader(abs_fits, 'H1_4861A')

    for j, line in enumerate(abs_target_lines):

        ion_array, wave_array, latex_array = lime.label_decomposition(line, scalar_output=True)

        cont_emis = fits.getdata(objFolder/'cont.fits', line)
        cont_abs = fits.getdata(objFolder/'abs_cont.fits', line)

        norm = cont_emis/cont_abs

        line_flux_emis = fits.getdata(emis_fits, line)
        line_flux_abs = fits.getdata(abs_fits, line) * -1

        # idxY, idxX = 167, 167
        # norm[idxY, idxX]
        # line_flux_emis[idxY, idxX]
        # line_flux_abs[idxY, idxX]
        # coeff_im[idxY, idxX]

        coeff_im = (line_flux_emis/(line_flux_abs * norm))

        print(np.sum(coeff_im<0))

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(projection=WCS(hdr_plot), slices=('x', 'y'))
        im = ax.imshow(coeff_im, norm=colors.LogNorm())
        cbar = fig.colorbar(im, ax=ax)
        title = f'{latex_array} ' + r'$\frac{emission}{absorption}$'
        ax.update({'title': r'Galaxy {}: {}'.format(obj, title), 'xlabel': r'RA', 'ylabel': r'DEC'})
        ax.set_xlim(95, 205)
        ax.set_ylim(75, 225)
        plt.show()
        # plt.savefig(objFolder/'line_ratios'/f'map_{obj}_OIII_ratio.png')