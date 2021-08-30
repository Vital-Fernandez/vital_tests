import numpy as np
import pandas as pd
from pathlib import Path
import src.specsiser as sr

from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, patches, _pylab_helpers
from astropy.wcs import WCS


def change_matplotlib_font(fig_input, prop_legend):
    for ax_obj in fig_input.canvas.figure.get_axes():
        ax_obj.legend(prop=prop_legend)
        print('zero')
        for label in ax_obj.get_xticklabels():
            print('un')
            label.set_fontproperties(prop_legend)
        for label in ax_obj.get_yticklabels():
            print('dos')
            label.set_fontproperties(prop_legend)

# Set figure conf

c_KMpS = 299792.458  # Speed of light in Km/s (https://en.wikipedia.org/wiki/Speed_of_light)

data_folder = Path('/home/vital/Astro-data/Observations/J0925 - OIII')

# file_list = ['J0925_OIII_ISIS.fits', 'J0925_OIII_MEGARA.fits', 'J0925_OIII_XSHOOTER_IZOTOV_2.fits', 'J0925_OIII_XSHOOTER_PUBLIC.fits']
file_list = ['J0925_OIII_ISIS.fits', 'J0925_OIII_MEGARA.fits', 'J0925_OIII_XSHOOTER_IZOTOV_2.fits']

spec_dict = {}
fit_dict = {}

z_obj = 0.3012097
waves_array = np.array([6480, 6495, 6500, 6530, 6540, 6548])
waves_array = waves_array / (1 + z_obj)

# fit_conf = {'O3_5007A_b': 'O3_5007A-O3_5007A_w1-O3_5007A_w2',
#             'O3_5007A_w1_sigma': {'expr': '>1.0*O3_5007A_sigma'},
#             'O3_5007A_w2_sigma': {'expr': '>1.0*O3_5007A_w1_sigma'},
#             'O3_5007A_cont_slope': {'vary': False},
#             'O3_5007A_cont_intercept': {'vary': False}}

fit_conf = {'O3_5007A_b': 'O3_5007A-O3_5007A_w1',
            'O3_5007A_w1_sigma': {'expr': '>1.0*O3_5007A_sigma'},
            'O3_5007A_cont_slope': {'vary': False},
            'O3_5007A_cont_intercept': {'vary': False}}


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
    # lm.plot_line_velocity(output_address=data_folder/f'{instrument}_velocity_plot', dark_mode=False)
    # plt.show()
    # lm.save_lineslog(lm.linesDF, data_folder/f'{instrument}_measurements_log.txt')

    print(f'Treating: {file} with {instrument} configuration')
    spec_dict[instrument] = [wave, flux]
    fit_dict[file] = lm


sizing_dict = {}
sizing_dict['font.family'] = 'Times New Roman'
sizing_dict['figure.figsize'] = (12, 10)
sizing_dict['axes.labelsize'] = 20
sizing_dict['axes.titlesize'] = 20
sizing_dict['xtick.labelsize'] = 16
sizing_dict['ytick.labelsize'] = 16
sizing_dict['legend.fontsize'] = 16


rcParams.update(sizing_dict)

lineLabel = 'O3_5007A'
fig, ax = plt.subplots()
for instr, lm_obj in fit_dict.items():

    idcsEmis, idcsContBlue, idcsContRed = lm_obj.define_masks(lm_obj.wave_rest,
                                                                lm_obj.flux,
                                                                lm_obj.lineWaves,
                                                                merge_continua=False)

    vel_plot = c_KMpS * (lm_obj.wave[idcsEmis] - lm_obj.peak_wave) / lm_obj.peak_wave
    cont_plot = lm_obj.m_cont * lm_obj.wave[idcsEmis] + lm_obj.n_cont
    flux_plot = lm_obj.flux[idcsEmis] - cont_plot

    ax.plot(vel_plot, flux_plot/np.max(flux_plot), label=instr)

    if 'ISIS' in instr:

        trans = ax.get_xaxis_transform()

        for i_percentil, percentil in enumerate([5, 10, 50, 90, 95]):
            vel_percentil = lm_obj.linesDF.loc[lineLabel, f'v_{percentil}']
            label_text = None if i_percentil > 0 else r'$v_{P_{i}}$'
            label_plot = r'$v_{{{}}}$'.format(percentil)
            ax.axvline(x=vel_percentil, label=label_text, linestyle='dotted', alpha=0.5, color='black')
            ax.text(vel_percentil, 0.650, label_plot, ha='center', va='center',
                    rotation='vertical', transform=trans)

        v_peak = 0

        v_med = lm_obj.linesDF.loc[lineLabel, f'v_med']
        label_vmed = r'$v_{{med}}={:0.1f}\,km/s$'.format(v_med)
        ax.axvline(x=v_med, label=label_vmed, color='black', linestyle='dashed', alpha=0.5)

        FWHM_int = lm_obj.linesDF.loc[lineLabel, f'FWHM_int']
        label_arrow = r'$FWHM={:0.1f}\,km/s$'.format(FWHM_int)
        p2 = patches.FancyArrowPatch((-FWHM_int/2, 0.5),
                                     (FWHM_int/2, 0.5),
                                     label=label_arrow,
                                     arrowstyle='<->',
                                     color='tab:red',
                                     transform=trans,
                                     mutation_scale=20)
        ax.add_patch(p2)

        v_5, v_95 = lm_obj.linesDF.loc[lineLabel, f'v_5'], lm_obj.linesDF.loc[lineLabel, f'v_95']
        v_10, v_90 = lm_obj.linesDF.loc[lineLabel, f'v_10'], lm_obj.linesDF.loc[lineLabel, f'v_90']
        w80 = v_90 - v_10
        w90 = v_95 - v_5

        label_arrow = r'$w_{{80}}={:0.1f}\,km/s$'.format(w80)
        p1 = patches.FancyArrowPatch((v_10, 0.4),
                                     (v_90, 0.4),
                                     label=label_arrow,
                                     arrowstyle='<->',
                                     color='tab:blue',
                                     transform=trans,
                                     mutation_scale=20)
        ax.add_patch(p1)

        A = ((v_90 - v_med) - (v_med-v_10)) / w80
        K = w90 / (1.397 * FWHM_int)

        label_A = r'$A = {:0.1f}$'.format(A)
        label_K = r'$K = {:0.1f}$'.format(K)
        ax.scatter(0, 1, label=label_A, alpha=0)
        ax.scatter(0, 1, label=label_K, alpha=0)

    # print(lm_obj.linesDF.loc['O3_5007A'].intg_flux)


ax.legend(prop={"family": "Times New Roman"})
# change_matplotlib_font(fig, prop_legend={"family": "Times New Roman"})
ax.update({'xlabel': 'Velocity (km/s)', 'ylabel': 'Normalized flux'})
plt.tight_layout()
plt.show()

# fig, ax = plt.subplots()
# for intr, data in spec_dict.items():
#     wave_plot, flux_plot = data
#     ax.plot(wave_plot, flux_plot/np.max(flux_plot), label=intr)
# ax.legend()
# ax.update({'xlabel': 'Wavelength', 'ylabel': 'Flux normalized by peak value', 'title': 'Instrument comparison'})
# plt.show()

