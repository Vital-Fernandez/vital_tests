import numpy as np
import src.specsiser as sr
import time
import lime

from pathlib import Path
from astro.papers.muse_CGCG007.muse_CGCG007_methods import chemical_lines_indexing, import_muse_fits
from astropy.io import fits
from astropy.table import Table
from lime.plots import STANDARD_PLOT
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from src.specsiser.plots import latex_labels
import pyneb as pn


a = np.array([[1, 2, 3],
              [4, np.nan, 6],
              [7, 8, 9]])

b = np.array([[1, 2, 3],
              [np.nan, 5, np.nan],
              [7, 8, 9]])

idcs_a = np.isnan(a)
idcs_b = np.isnan(b)

print(idcs_a | idcs_b)


# idcs_nan = ~np.isnan(a)
# vector = a[idcs_nan].flatten()
#
# print(vector)
#
# a[idcs_nan] = vector * 2
#
# print(a)

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

merge_dict = {'O2_7319A_b': 'O2_7319A-O2_7330A'}

R_v = obsData['Extinction']['R_v']
red_law = obsData['Extinction']['red_law']
image_size = obsData['sample_data']['grid_shape_array'].astype(int)

lines_mask_dict = {'MASK_0': ['H1_4861A', 'H1_9229A', 'H1_9015A', 'H1_8863A', 'H1_8750A'],
                   'MASK_1': ['H1_4861A', 'H1_9229A', 'H1_9015A', 'H1_8863A'],
                   'MASK_2': ['H1_4861A', 'H1_6563A', 'H1_9229A', 'H1_9015A'],
                   'MASK_3': ['H1_4861A', 'H1_6563A'],
                   'MASK_4': ['H1_4861A', 'H1_6563A'],
                   'MASK_5': ['H1_4861A', 'H1_6563A']}

S2 = pn.Atom('S', 2)

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    chemFolder = f'{objFolder}/chemistry'
    obsLog_addresss = objFolder / f'{obj}_linesLog.fits'
    absLog_addresss = objFolder / f'{obj}_linesLog_abs.fits'
    maskFits_address = objFolder/f'{obj}_masks.fits'
    db_address = objFolder / f'{obj}_database.fits'

    # Output files
    output_map = objFolder/f'{obj}_neSII_gauss.fits'
    images_dict = {'nSII': np.full(image_size, np.nan),
                   'nSII_err': np.full(image_size, np.nan)}

    # Header data with the cosmological coordinates
    hdr_cube = fits.getheader(db_address, extname='PlotConf')
    hdr_plot = {}
    for entry in lime.COORD_ENTRIES:
        if entry in hdr_cube:
            hdr_plot[entry] = hdr_cube[entry]

    S2_6716A_flux, S2_6716A_err = fits.getdata(objFolder/'gauss_flux.fits', extname='S2_6716A'), fits.getdata(objFolder/'gauss_err.fits', extname='S2_6716A')
    S2_6731A_flux, S2_6731A_err = fits.getdata(objFolder/'gauss_flux.fits', extname='S2_6731A'), fits.getdata(objFolder/'gauss_err.fits', extname='S2_6731A')

    idcs_nan = np.isnan(S2_6716A_flux) | np.isnan(S2_6731A_flux)
    SII_ratio = S2_6716A_flux[~idcs_nan]/S2_6731A_flux[~idcs_nan]

    print(f'TODOS NAN nSII_vector: {np.all(np.isnan(SII_ratio.flatten()))}')

    nSII_vector = S2.getTemDen(SII_ratio.flatten(), tem=10000.0, to_eval='L(6716)/L(6731)')
    nSII_map = np.full(S2_6716A_flux.shape, np.nan)
    nSII_map[~idcs_nan] = nSII_vector

    print(f'TODOS NAN nSII_vector: {np.all(np.isnan(nSII_vector))}')
    print(f'TODOS NAN nSII_map: {np.all(np.isnan(nSII_map))}')

    paramHDUs = fits.HDUList()
    paramHDUs.append(fits.PrimaryHDU())
    paramHDUs.append(fits.ImageHDU(name='ne_SII', data=nSII_map, header=fits.Header(hdr_plot), ver=1))
    paramHDUs.writeto(output_map, overwrite=True, output_verify='fix')


    # Write to new file




    #     # Primary header
    #     paramHDUs = fits.HDUList()
    #     paramHDUs.append(fits.PrimaryHDU())
    #
    #     # For param in image for the parameter maps
    #     for param, image in images_dict.items():
    #         hdr = fits.Header({'PARAM': param})
    #         hdr.update(plot_dict)
    #         paramHDUs.append(fits.ImageHDU(name=param, data=image, header=hdr, ver=1))
    #
    #     # Write to new file
    #     paramHDUs.writeto(output_map, overwrite=True, output_verify='fix')
    #


    # plot_dict = {'CRPIX1': hdr_plot['CRPIX1'],
    #              'CRPIX2': hdr_plot['CRPIX2'],
    #              'CD1_1': hdr_plot['CD1_1'],
    #              'CD1_2': hdr_plot['CD1_2'],
    #              'CD2_1': hdr_plot['CD2_1'],
    #              'CD2_2': hdr_plot['CD2_2'],
    #              'CUNIT1': hdr_plot['CUNIT1'],
    #              'CUNIT2': hdr_plot['CUNIT2'],
    #              'CTYPE1': hdr_plot['CTYPE1'],
    #              'CTYPE2': hdr_plot['CTYPE2']}
    #
    # # Open the fits
    # obsHDUs = fits.open(obsLog_addresss, lazy_load_hdus=False)
    # absHDUs = fits.open(absLog_addresss, lazy_load_hdus=False)
    #
    # # Loop throught the line regions
    # start = time.time()
    # for idx_region in [0, 1, 2, 3, 4, 5]:
    #
    #     # Voxel mask
    #     region_label = f'mask_{idx_region}'
    #     region_mask = fits.getdata(maskFits_address, region_label, ver=1)
    #     region_mask = region_mask.astype(bool)
    #     idcs_voxels = np.argwhere(region_mask)
    #
    #     # Simulation time estimation
    #     n_voxels = idcs_voxels.shape[0]
    #     print(f'\nTreating {region_label} consisting of {n_voxels}')
    #
    #     # Region chemical configuration
    #     # chem_conf_file = dataFolder/f'{obj}_chemical_model_region_{idx_region}.txt'
    #     # chem_conf = lime.load_cfg(chem_conf_file)
    #     # chem_conf = lime.load_cfg(chem_conf_file)
    #
    #     # Extinction law
    #     red = sr.ExtinctionModel(Rv=R_v, red_curve=red_law)
    #     input_lines = lines_mask_dict[f'MASK_{idx_region}']
    #
    #     for idx_voxel, idx_pair in enumerate(idcs_voxels):
    #         idx_j, idx_i = idx_pair
    #         print(f'\nTreating voxel {idx_j}-{idx_i}: ({idx_voxel}/{n_voxels})')
    #
    #         # Load voxel data:
    #         ext_ref = f'{idx_j}-{idx_i}_linelog'
    #         # obs_log = lime.load_lines_log(obsLog_addresss, ext_ref)
    #         # abs_log = lime.load_lines_log(absLog_addresss, ext_ref)
    #
    #         if ext_ref in obsHDUs:
    #             obs_log = Table.read(obsHDUs[ext_ref]).to_pandas()
    #             obs_log.set_index('index', inplace=True)
    #         else:
    #             obs_log = None
    #
    #         if ext_ref in absHDUs:
    #             abs_log = Table.read(absHDUs[ext_ref]).to_pandas()
    #             abs_log.set_index('index', inplace=True)
    #         else:
    #             abs_log = None
    #
    #         if obs_log is not None and abs_log is not None:
    #
    #             if 'H1_4861A' in obs_log.index:
    #
    #                 # Establish and normalize the lines we want
    #                 linesDF = chemical_lines_indexing(input_lines, obs_log, abs_log, obsData, recomb_all=False)
    #
    #                 if len(linesDF.index) > 1:
    #                     cHbeta, cHbeta_err = red.cHbeta_from_log(linesDF)#, plot_address=True)
    #                     if np.isnan(cHbeta):
    #                         print('ONLY ONE LINE')
    #
    #                     images_dict['cHbeta'][idx_j, idx_i] = cHbeta
    #                     images_dict['cHbeta_err'][idx_j, idx_i] = cHbeta_err
    #
    #     # Primary header
    #     paramHDUs = fits.HDUList()
    #     paramHDUs.append(fits.PrimaryHDU())
    #
    #     # For param in image for the parameter maps
    #     for param, image in images_dict.items():
    #         hdr = fits.Header({'PARAM': param})
    #         hdr.update(plot_dict)
    #         paramHDUs.append(fits.ImageHDU(name=param, data=image, header=hdr, ver=1))
    #
    #     # Write to new file
    #     paramHDUs.writeto(output_map, overwrite=True, output_verify='fix')
    #
    # obsHDUs.close()
    # absHDUs.close()

    # # ----------------------------------------- Generate the image plots ----------------------------------------
    # flux6563_image = fits.getdata(db_address, f'H1_6563A_flux', ver=1)
    # halpha_min_level = fits.getval(db_address, keyword=f'P9050', extname=f'H1_6563A_flux')
    # halpha_thresd_level = fits.getval(db_address, keyword=f'P9250', extname=f'H1_6563A_flux')
    #
    # with fits.open(output_map) as hdu_list:
    #     param = 'cHbeta'
    #     image_data, image_header = hdu_list[param].data, hdu_list[param].header
    #
    #     # idcx_high = image_data > 3.0
    #     # image_data[idcx_high] = np.nan
    #     print(np.sum(~np.isnan(image_data)))
    #
    #     defaultConf = STANDARD_PLOT.copy()
    #     rcParams.update(defaultConf)
    #
    #     halpha_cmap = cm.gray
    #     halpha_cmap.set_under('black')
    #
    #     fig = plt.figure(figsize=(10, 10))
    #     ax = fig.add_subplot(projection=WCS(image_header), slices=('x', 'y'))
    #
    #     bg_color = colors.SymLogNorm(linthresh=halpha_thresd_level, vmin=halpha_min_level, base=10)
    #     im = ax.imshow(flux6563_image, cmap=halpha_cmap, norm=bg_color)
    #
    #     param_min, param_max = np.nanmin(image_data), np.nanmax(image_data)
    #     divnorm = colors.TwoSlopeNorm(vcenter=0.0, vmin=-0.05, vmax=0.5, )
    #
    #     # im2 = ax.imshow(image_data, norm=colors.LogNorm())
    #     im2 = ax.imshow(image_data, norm=divnorm)
    #
    #     cbar = fig.colorbar(im2, ax=ax)
    #     param_label = latex_labels[param]
    #     ax.update({'title': r'Galaxy {}, {}'.format(obj, param_label), 'xlabel': r'RA', 'ylabel': r'DEC'})
    #     ax.set_xlim(95, 210)
    #     ax.set_ylim(80, 222)
    #     # plt.savefig(objFolder/f'{obj}_parameter_map_{ext_label}')
    #     plt.show()


