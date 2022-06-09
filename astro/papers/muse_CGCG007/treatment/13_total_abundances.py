import numpy as np
import lime
from lime.plots import STANDARD_PLOT
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import save_log_maps, latex_labels, signif_figures, param_units
from lime.io import save_cfg


def to_natural_abund(abund_array):
    return np.power(10, abund_array - 12)


def to_log_abund(abund_array):
    return 12 + np.log10(abund_array)


# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']

fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

voxel_grid_size = obsData['sample_data']['grid_shape_array']

for i, obj in enumerate(objList):

    # Data location
    objFolder = resultsFolder/obj
    maskFits_address, mask_list = objFolder/f'{obj}_masks.fits', ['MASK_0', 'MASK_1', 'MASK_2']
    db_address = objFolder/f'{obj}_database.fits'
    outputDb = objFolder/f'{obj}_chemical.fits'
    chemFolder = objFolder/'chemistry'

    # Parameters to plot
    param_list = np.array(['n_e', 'T_low', 'cHbeta', 'Ar4', 'Ar3', 'O2', 'O3', 'N2', 'He1', 'S2', 'S3'])

    # Get mask indeces:
    spatial_mask_dict = {}
    with fits.open(maskFits_address) as hdu_masks:
        for mask_name in mask_list:
            mask_data = hdu_masks[mask_name].data.astype(bool)
            spatial_mask_dict[mask_name] = mask_data
    total_mask = np.array(list(spatial_mask_dict.values()))
    total_mask = total_mask.sum(axis=0).astype(bool)

    # # Data for the astronomical coordinates
    # hdr = fits.getheader(maskFits_address, extname='MASK_0')
    #
    # # Generate the map files
    # save_log_maps(outputDb, param_list, chemFolder, maskFits_address, mask_list, ext_log='_CHEMISTRY_OUTPUTS',
    #               page_hdr=hdr)

    # ----------------------------------------- Generate the parameter histograms ----------------------------------------
    store_dict, chain_dict = {}, {}
    for param in param_list:

        with fits.open(f'{chemFolder}/{param}.fits') as hdu_list:

            image_data, image_header = hdu_list[param].data, hdu_list[param].header
            err_data = hdu_list[f'{param}_err'].data

            array_data = image_data[total_mask]
            array_err = err_data[total_mask]

            chain_dict[param] = array_data
            print(f'{param}: {np.nanmean(chain_dict[param])}')

    # Oxygen abundance
    O2, O3 = to_natural_abund(chain_dict['O2']), to_natural_abund(chain_dict['O3'])
    OH = O2 + O3

    # Nirogen abundance
    N2 = to_natural_abund(chain_dict['N2'])
    NO = N2/O2
    NH = NO * OH

    # Argon abundance
    Ar3, Ar4 = to_natural_abund(chain_dict['Ar3']), to_natural_abund(chain_dict['Ar4'])
    ArH = Ar3 + Ar4

    # Sulfur abundance
    S2, S3 = to_natural_abund(chain_dict['S2']), to_natural_abund(chain_dict['S3'])
    m_conf, n_conf = np.random.normal(1.162, 0.006, Ar3.size), np.random.normal(0.05, 0.01, Ar3.size)
    exp_value = (np.log10(Ar3/Ar4) - n_conf) / m_conf
    S3S4 = np.power(10, exp_value)
    S4 = S3/S3S4
    SH = S2 + S3 + S4
    ICF_S4 = SH/(S2 + S3)
    SO = SH/OH

    #Extra params
    O2_O3 = O2/O3
    S2_S3 = S2/S3
    eta = O2_O3/S2_S3

    # Save mean values to log
    He1 = chain_dict['He1']
    SO_dors = np.power(10, -1.78), np.power(10, -0.02)
    SO_dist = np.random.normal(SO_dors[0], SO_dors[1], He1.size)
    Y_O = (4 * He1 * (1 - 20 * OH)) / (1 + 4*He1)
    Y_S = (4 * He1 * (1 - 20 * SO_dist * SH)) / (1 + 4*He1)

    # Store to a dictiory
    param_chain_dict = dict(OH=to_log_abund(OH),
                      NO=np.log10(NO),
                      NH=to_log_abund(NH),
                      ArH=to_log_abund(ArH),
                      S4=to_log_abund(S4),
                      SH=to_log_abund(SH),
                      ICF_S4=ICF_S4,
                      SO=np.log10(SO),
                      Y_O=Y_O,
                      Y_S=Y_S,
                      S2_S3=S2_S3,
                      O2_O3=O2_O3,
                      eta=eta)

    # Save total abundance distributions to the configuration file
    store_dict = {}
    print(f'Esto es {np.nanmedian(param_chain_dict["OH"])}')
    for param, chain in param_chain_dict.items():
        median = np.round(np.nanmedian(chain), signif_figures[param])
        upper_limit = np.round(np.nanpercentile(chain, 84) - np.nanmedian(chain), signif_figures[param]+1)
        lower_limit = np.round(np.nanmedian(chain) - np.nanpercentile(chain, 16), signif_figures[param]+1)
        n_voxels = np.sum(~np.isnan(chain))
        store_dict[f'{param}_array'] = np.array([median, upper_limit, lower_limit, n_voxels])
        print(param, median, upper_limit, lower_limit, n_voxels)
    save_cfg('../muse_CGCG007.ini', store_dict, section_name='Total_abundances_direct_method', clear_section=True)
