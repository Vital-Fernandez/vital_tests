import numpy as np
import pandas as pd
from pathlib import Path
import src.specsiser as sr
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams
from astropy.wcs import WCS

# Set figure conf
sizing_dict = {}
sizing_dict['figure.figsize'] = (20, 8)
sizing_dict['legend.fontsize'] = 12
sizing_dict['axes.labelsize'] = 12
sizing_dict['axes.titlesize'] = 12
sizing_dict['xtick.labelsize'] = 12
sizing_dict['ytick.labelsize'] = 12

rcParams.update(sizing_dict)

data_folder = Path('/home/vital/Astro-data/Observations/J0925 - OIII')

file_list = ['J0925_OIII_ISIS.fits', 'J0925_OIII_MEGARA.fits', 'J0925_OIII_XSHOOTER_IZOTOV_2.fits', 'J0925_OIII_XSHOOTER_PUBLIC.fits']

data_dict = {}

z_obj = 0.3012097
waves_array = np.array([6480, 6495, 6500, 6530, 6540, 6548])
waves_array = waves_array / (1 + z_obj)

fit_conf = {'O3_5007A_b': 'O3_5007A-O3_5007A_w1-O3_5007A_w2',
            'O3_5007A_w1_sigma': {'expr': '>1.0*O3_5007A_sigma'},
            'O3_5007A_w2_sigma': {'expr': '>1.0*O3_5007A_w1_sigma'},
            'O3_5007A_cont_slope': {'vary': False},
            'O3_5007A_cont_intercept': {'vary': False}}

# fit_conf = {'O3_5007A_b': 'O3_5007A-O3_5007A_w1',
#             'O3_5007A_w1_sigma': {'expr': '>1.0*O3_5007A_sigma'},
#             'O3_5007A_cont_slope': {'vary': False},
#             'O3_5007A_cont_intercept': {'vary': False}}

z_obj = 0.3012097
for i, file in enumerate(file_list):

    if 'ISIS' in file:
        instrument, ext = 'ISIS', 0

        with fits.open(data_folder/file) as hdul:
            data, header = hdul[0].data, hdul[0].header
            wcs4 = WCS(hdul[0].header)

        index4 = np.arange(header['NAXIS1'])
        wave = wcs4.wcs_pix2world(index4, 0, 0, 0)[0]

        w_min = header['CRVAL1']
        dw = header['CD1_1']  # dw = 0.862936 INDEF (Wavelength interval per pixel)
        pixels = header['NAXIS1']  # nw = 3801 number of output pixels
        w_max = w_min + dw * pixels
        # wave = np.linspace(w_min, w_max, pixels, endpoint=False)
        flux = data

    if 'XSHOOTER_PUBLIC' in file:
        instrument, ext = 'xshooter-public', 1

        with fits.open(data_folder/file) as hdul:
            data, header = hdul[1].data, hdul[1].header

        wave = data['WAVE'][0] * 10
        flux = data['FLUX'][0]
        err = data['ERR'][0]
        crop_array = np.array([6470, 6550])
        idcs_crop = np.searchsorted(wave, crop_array)
        wave, flux, err = wave[idcs_crop[0]:idcs_crop[1]], flux[idcs_crop[0]:idcs_crop[1]], err[idcs_crop[0]:idcs_crop[1]]

    if 'XSHOOTER_IZOTOV' in file:

        instrument, ext = 'XSHOOTER-Izotov', 1

        with fits.open(data_folder/file) as hdul:
            data, header = hdul[0].data, hdul[0].header
            wcs4 = WCS(header)

        index4 = np.arange(header['NAXIS1'])
        wave = wcs4.wcs_pix2world(index4, 0)[0]
        flux = data

        crop_array = np.array([6470, 6550])
        idcs_crop = np.searchsorted(wave, crop_array)
        wave, flux = wave[idcs_crop[0]:idcs_crop[1]], flux[idcs_crop[0]:idcs_crop[1]]


    if 'MEGARA' in file:

        instrument = 'MEGARA'

        with fits.open(data_folder/file) as hdul:
            data, header = hdul[0].data, hdul[0].header
            wcs4 = WCS(header)

        index4 = np.arange(header['NAXIS1'])
        wave = wcs4.wcs_pix2world(index4, 0)[0]

        w_min = header['CRVAL1']
        dw = header['CDELT1']  # dw (Wavelength interval per pixel)
        pixels = header['NAXIS1']  # nw number of output pixels
        w_max = w_min + dw * pixels
        # wave = np.linspace(w_min, w_max, pixels, endpoint=False)
        flux = data

    lm = sr.LineMesurer(wave, flux, redshift=z_obj, normFlux=np.median(flux))
    lm.fit_from_wavelengths('O3_5007A_b', line_wavelengths=waves_array, user_conf=fit_conf)
    # lm.print_results(show_plot=True)
    lm.plot_line_velocity(output_address=data_folder/f'{instrument}_velocity_plot', dark_mode=False)
    # lm.save_lineslog(lm.linesDF, data_folder/f'{instrument}_measurements_log.txt')

    print(f'Treating: {file} with {instrument} configuration')
    data_dict[instrument] = [wave, flux/np.max(flux)]


fig, ax = plt.subplots()
for intr, data in data_dict.items():
    wave_plot, flux_plot = data
    ax.plot(wave_plot, flux_plot, label=intr)
ax.legend()
ax.update({'xlabel': 'Wavelength', 'ylabel': 'Flux normalized by peak value', 'title': 'Instrument comparison'})
plt.show()

6515.157