import numpy as np
import lime
from lime.plots import STANDARD_PLOT
from pathlib import Path
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from astro.papers.muse_CGCG007.muse_CGCG007_methods import save_log_maps, latex_labels, signif_figures, param_units
from lime.io import save_cfg
from bces_code import bcesp
from uncertainties import unumpy, ufloat
from scipy.optimize import curve_fit
from scipy.odr import   Model, ODR, RealData

def to_natural_abund(abund_array):
    return np.power(10, abund_array - 12)


def to_log_abund(abund_array):
    return 12 + np.log10(abund_array)


def linear_model(x, m, n):
    return m * x + n

def f(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]

# Declare data and files location
obsData = lime.load_cfg('../muse_CGCG007.ini')
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']

# fitsFolder = Path(obsData['data_location']['fits_folder'])
# dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path('D:/Dropbox/Astrophysics/Papers/muse_CGCG007/treatment')

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

    # ----------------------------------------- Generate the parameter histograms ----------------------------------------
    store_dict, chain_dict = {}, {}
    for param in param_list:

        with fits.open(f'{chemFolder}/{param}.fits') as hdu_list:

            image_data, image_header = hdu_list[param].data, hdu_list[param].header
            err_data = hdu_list[f'{param}_err'].data

            array_data = image_data[total_mask]
            array_err = err_data[total_mask]

            chain_dict[param] = array_data
            chain_dict[f'{param}_err'] = array_err

            print(f'{param}: {np.nanmean(chain_dict[param])}')

    # Oxygen abundance
    He1_nom, He1_err = chain_dict['He1'], chain_dict['He1_err']
    O2_log, O2_err = chain_dict['O2'], chain_dict['O2_err']
    O3_log, O3_err = chain_dict['O3'], chain_dict['O3_err']
    O2 = unumpy.uarray(O2_log, O2_err)
    O3 = unumpy.uarray(O3_log, O3_err)
    He1 = unumpy.uarray(He1_nom, He1_err)

    OH = to_natural_abund(O2) + to_natural_abund(O3)

    # Save mean values to log
    SO_dors = ufloat(np.power(10, -1.78), np.power(10, -0.02))
    Y = (4 * He1 * (1 - 20 * OH)) / (1 + 4*He1)

    OH_nom, OH_err = unumpy.nominal_values(OH), unumpy.std_devs(OH)
    Y_nom, Y_err = unumpy.nominal_values(Y), unumpy.std_devs(Y)

    idcs_nan = ~np.isnan(OH_nom) & ~np.isnan(Y_nom)# & (OH_nom < 0.0004)
    OH_nom, OH_err = OH_nom[idcs_nan], OH_err[idcs_nan]
    Y_nom, Y_err = Y_nom[idcs_nan], Y_err[idcs_nan]

    linear_model = Model(f)
    mydata = RealData(OH_nom, Y_nom, sx=OH_err, sy=Y_err)
    myodr = ODR(mydata, linear_model, beta0=[0, 24.0])
    out = myodr.run()
    out.pprint()

    fig, ax = plt.subplots()
    ax.scatter(OH_nom[idcs_nan], Y_nom[idcs_nan])
    ax.errorbar(OH_nom[idcs_nan], Y_nom[idcs_nan], xerr=OH_err[idcs_nan], yerr=Y_err[idcs_nan], ls='none')
    ax.legend()
    ax.update({'xlabel': r'$\frac{O}{H}$', 'ylabel': 'Y', 'title': 'Gaussian fitting'})
    plt.show()

    # m, n, m_err, n_err, cov = bcesp(unumpy.nominal_values(OH_nom), unumpy.std_devs(OH_err),
    #                                     unumpy.nominal_values(Y_nom), unumpy.std_devs(Y_err),
    #                                     cerr=np.zeros(len(OH_nom)),
    #                                     nsim=100)
    #
    # print('BCES y dependent ')
    # print('n', n[0], n_err[0])
    # print('m', m[0], m_err[0])
    # print('BCES Orthogonal')
    # print('n', n[3], n_err[3])
    # print('m', m[3], m_err[3])

    # Curvefit
    best_vals, covar = curve_fit(linear_model, OH_nom, Y_nom)
    print(best_vals)



