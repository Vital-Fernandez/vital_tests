from pathlib import Path
from astropy.io import fits
import lime
import numpy as np

def getSpectra(specIndex, mikeArm):
    waveArm = mikeArm[1].data
    fluxArm = mikeArm[2].data
    errorArm = mikeArm[3].data
    W = waveArm[specIndex]
    F = fluxArm[specIndex]
    E = errorArm[specIndex]
    return W, F, E

# file_address = Path('/home/vital/Dropbox/Astrophysics/Data/CGCG007_mike/CGCG007_redArm_v2.fits')
file_address = Path('/home/vital/Dropbox/Astrophysics/Data/CGCG007_mike/CGCG007_redArm_v3.fits')
# file_address = Path('/home/vital/Dropbox/Astrophysics/Data/CGCG007_mike/CGCG007_redArm_v4.fits')

ext = 15
with fits.open(file_address) as hdu_list:
    wave_m = hdu_list[1].data[ext]
    flux_m = hdu_list[2].data[ext]
    err_m = hdu_list[3].data[ext]

spec_order = lime.Spectrum(wave_m, flux_m, err_m, redshift=0.004691)
spec_order.plot_spectrum(frame='rest')

line = 'H1_6563A_b'
lineWaves = np.array([6539, 6545, 6548, 6587, 6589, 6597])
fit_conf = {'H1_6563A_b': 'H1_6563A-H1_6563A_w1-H1_6563A_w2-N2_6584A-N2_6548A',
            'H1_6563A_w1_sigma': {'expr': '>2*H1_6563A_sigma'},
            'H1_6563A_w2_sigma': {'expr': '>3*H1_6563A_w1_sigma'},
            'H1_6563A_cont_slope': {'vary': False},
            'H1_6563A_cont_intercept': {'vary': False}}



spec_order.fit_from_wavelengths(line, lineWaves, user_cfg=fit_conf)
spec_order.display_results(show_fit_report=True, show_plot=True, frame='rest')


# file_address = Path('/home/vital/Dropbox/Astrophysics/Data/CGCG007_mike/CGCG007_redArm_v2.fits')
#
# ext = 14
# with fits.open(file_address) as hdu_list:
#     wave_m = hdu_list[1].data[ext]
#     flux_m = hdu_list[2].data[ext]
#     err_m = hdu_list[3].data[ext]
#
# spec_order = lime.Spectrum(wave_m, flux_m, err_m, redshift=0.004691)
# spec_order.plot_spectrum(frame='rest')
#
# line = 'H1_6563A_b'
# # lineWaves = np.array([6551, 6555, 6556, 6572, 6575, 6583])
# lineWaves = np.array([6539, 6545, 6548, 6587, 6589, 6597])
# fit_conf = {'H1_6563A_b': 'H1_6563A-H1_6563A_w1-H1_6563A_w2-N2_6584A-N2_6548A',
#             'H1_6563A_w1_sigma': {'expr': '>2*H1_6563A_sigma'},
#             'H1_6563A_w2_sigma': {'expr': '>3*H1_6563A_w1_sigma'}}
#
# spec_order.fit_from_wavelengths(line, lineWaves, user_cfg=fit_conf)
# spec_order.display_results(show_fit_report=True, show_plot=True, frame='rest')