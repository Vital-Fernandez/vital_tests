import numpy as np

import lime
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors

from astro.papers.muse_CGCG007.muse_CGCG007_methods import abs_target_lines

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
voxel_grid_size = obsData['sample_data']['grid_shape_array']
coordinates_keys_list = obsData['data_location']['wcs_key_list']

# ------ Emission/absroption distribution
for i, obj in enumerate(objList):

    # Data location
    cube_address = fitsFolder/fileList[i]
    objFolder = resultsFolder/obj
    db_address = objFolder / f'{obj}_database.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'

    abs_fits = objFolder/'abs_gauss_flux.fits'
    Hbeta_abs = fits.getdata(abs_fits, 'H1_4861A') * -1

    ion_Hbeta, wave_Hbeta, latex_Hbeta = lime.label_decomposition('H1_4861A', scalar_output=True)
    cont_emisHbeta = fits.getdata(objFolder / 'cont.fits', 'H1_4861A')
    cont_absHbeta = fits.getdata(objFolder / 'abs_cont.fits', 'H1_4861A')
    hdr_plot = fits.getheader(abs_fits, 'H1_4861A')

    colorNorm = colors.Normalize(0, 6)
    cmap = cm.get_cmap(name=None)
    abs_dict = {}
    abs_count = {}

    for j, line in enumerate(abs_target_lines):

        ion_array, wave_array, latex_array = lime.label_decomposition(line, scalar_output=True)

        array_container = []
        data_labels = []
        colors = []
        n_count = 0

        line_flux_abs = fits.getdata(abs_fits, line) * -1
        cont_emis = fits.getdata(objFolder / 'cont.fits', line)
        cont_abs = fits.getdata(objFolder / 'abs_cont.fits', line)

        for idx_region in [0, 1, 2, 3, 4, 5]:

            # Voxel mask
            region_label = f'MASK_{idx_region}'
            region_mask = fits.getdata(maskFits_address, region_label, ver=1)
            region_mask = region_mask.astype(bool)

            abs_norm = line_flux_abs[region_mask] * (cont_emis[region_mask]/cont_abs[region_mask])
            abs_Hbeta_norm = Hbeta_abs[region_mask] * (cont_emisHbeta[region_mask]/cont_absHbeta[region_mask])

            coeff_im = abs_norm/abs_Hbeta_norm
            idcs_non_nan = ~np.isnan(coeff_im) & (coeff_im > 0) & (coeff_im < 4)

            if idcs_non_nan.size > 0:
                data_array = coeff_im[idcs_non_nan]
                array_container.append(data_array)
                data_labels.append(f'{region_label} ({len(data_array)})')
                colors.append(cmap(colorNorm(idx_region)))
                n_count += len(data_array)

        # Computing mean absorptions
        data_total = np.concatenate(array_container)
        mean_abs, std_abs = data_total.mean(), data_total.std()
        abs_dict[f'{line}_abs_array'] = np.array([mean_abs, std_abs])
        abs_count[f'{line}_count'] = n_count

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.hist(array_container, label=data_labels, bins=30, log=True, stacked=True, color=colors)
        ratio_label = r'$\frac{{{}}}{{{}}}$'.format(latex_array.replace('$', ''), latex_Hbeta.replace('$', ''))
        title = f'{ratio_label} absorption'
        ax.axvline(x=mean_abs, color='black', linestyle='--')
        ax.axvspan(mean_abs - std_abs, mean_abs + std_abs, alpha=0.25, color='grey')
        ax.update({'title': r'Galaxy {}'.format(obj), 'xlabel': title, 'ylabel': r'Voxel count'})
        ax.legend()
        plt.show()
        # plt.savefig(f'{objFolder}/absorptions/{line}_abs.png')

    lime.io.save_cfg('../muse_CGCG007.ini', abs_dict, f'CGCG007_absorptions', clear_section=True)
    lime.io.save_cfg('../muse_CGCG007.ini', abs_count, f'CGCG007_abs_count', clear_section=True)


