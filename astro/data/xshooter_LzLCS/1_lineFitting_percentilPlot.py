import numpy as np
from astropy.io import fits
from pathlib import Path
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams, transforms, patches
from scipy.interpolate import interp1d

c_KMpS = 299792.458  # Speed of light in Km/s (https://en.wikipedia.org/wiki/Speed_of_light)


obsData = sr.loadConfData('./xshooter_LzLCS.ini')
data_folder = Path(obsData['data_location']['data_folder'])
results_folder = Path(obsData['data_location']['results_folder'])
objfile_list = obsData['data_location']['objfile_list']
sigmafile_list = obsData['data_location']['sigmafile_list']
objRef_list = obsData['data_location']['ref_list']
maskfile = obsData['data_location']['generalMask']

wmin_array = obsData['sample_data']['w_min_array']
wmax_array = obsData['sample_data']['w_max_array']
norm_flux = obsData['sample_data']['norm_flux']
z_obj = obsData['sample_data']['z_obj']
profile_conf = obsData['line_fitting']

verbose = True

for i, objName in enumerate(objRef_list):

    # Input data
    spec_file, sigm_file = data_folder/objfile_list[i], data_folder/sigmafile_list[i]

    # Output data
    lineslog_file = results_folder/f'{objName}_linesLog.txt'
    lineslog_table = results_folder/f'{objName}_flux_table'

    # Load inputs
    wave, flux, header = sr.import_fits_data(spec_file, instrument='xshooter', frame_idx=0)
    wave_sigma, sigma, header_sigma = sr.import_fits_data(sigm_file, instrument='xshooter', frame_idx=0)

    lm = sr.LineMesurer(wave, flux, crop_waves=[wmin_array[i], wmax_array[i]], input_err=sigma, normFlux=norm_flux, redshift=z_obj)

    # lm.plot_spectrum(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF, specLabel=f'{objName}')
    mask_local = data_folder/f'{objName}_mask'
    mask_local_df = sr.lineslogFile_to_DF(mask_local)

    obsLines = mask_local_df.index.values
    for j, lineLabel in enumerate(obsLines):
        wave_regions = mask_local_df.loc[lineLabel, 'w1':'w6'].values
        if lineLabel in ['O3_5007A_b']:
            # wave_regions[2] = 8400.0/(1+lm.redshift)
            # wave_regions[3] = 8445.0/(1+lm.redshift)
            lm.fit_from_wavelengths(lineLabel, wave_regions, user_conf=profile_conf)
            lm.print_results(show_plot=True, show_fit_report=False, log_scale=False, frame='rest')

            lm.plot_line_velocity()

            lm.linesDF
            # idcsEmis, idcsContBlue, idcsContRed = lm.define_masks(lm.wave_rest,
            #                                                       lm.flux,
            #                                                       wave_regions,
            #                                                       merge_continua=False)
            #
            # # x_plot = lm.wave[idcsEmis]
            # x_plot = c_KMpS * (lm.wave[idcsEmis]-lm.peak_wave)/lm.peak_wave
            # x_in = lm.wave[idcsEmis]
            # pix_size = np.diff(x_in).mean()
            # vel_plot = lm.flux[idcsEmis]
            # cont_plot = lm.m_cont * lm.wave[idcsEmis] + lm.n_cont
            #
            # percentile_array = np.zeros(vel_plot.size)
            # for i_pix in np.arange(vel_plot.size):
            #     i_flux = (vel_plot[:i_pix].sum() - cont_plot[:i_pix].sum()) * pix_size
            #     percentile_array[i_pix] = i_flux/lm.intg_flux * 100
            #
            # vel_med = np.median(vel_plot)
            #
            # target_percentiles = np.array([2, 5, 10, 50, 90, 95, 98])
            # Interpolation = interp1d(percentile_array, x_plot, kind='slinear')
            # vel_percentiles = Interpolation(target_percentiles)
            # fig, ax = plt.subplots()
            # ax.step(x_plot, vel_plot, label=r'{}'.format(lineLabel), where='mid')
            # trans = ax.get_xaxis_transform()
            #
            # for i_percentil, percentil in enumerate(target_percentiles):
            #     label_plot = r'$v_{{{}}}$'.format(percentil)
            #     label_text = None if i_percentil > 0 else r'$v_{Pth}$'
            #     ax.axvline(x=vel_percentiles[i_percentil], label=label_text, color='grey', linestyle='dotted', alpha=0.5)
            #     ax.text(vel_percentiles[i_percentil], 0.80, label_plot, ha='center', va='center',
            #             rotation='vertical', backgroundcolor='white', transform=trans, alpha=0.5)
            #
            # w80 = vel_percentiles[4]-vel_percentiles[2]
            # label_arrow = r'$w_{{80}}={:0.1f}\,Km/s$'.format(w80)
            # p1 = patches.FancyArrowPatch((vel_percentiles[2], 0.5),
            #                              (vel_percentiles[4], 0.5),
            #                              label=label_arrow,
            #                              arrowstyle='<->',
            #                              color='tab:blue',
            #                              transform=trans,
            #                              mutation_scale=20)
            # ax.add_patch(p1)
            # label_vmed = r'$v_{{med}}={:0.1f}\,Km/s$'.format(vel_med)
            # ax.axvline(x=vel_med, color='black', label=label_vmed, linestyle='dashed', alpha=0.5)
            #
            # label_vmed = r'$v_{{peak}}$'
            # ax.axvline(x=0.0, color='black', label=label_vmed, alpha=0.5)
            #
            #
            #
            # # ax.annotate(vel_percentiles[2], 0.5,
            # #          vel_percentiles[4]-vel_percentiles[2], 0,
            # #          color='red', head_length=0.07, head_width=0.025,
            # #          transform=trans,
            # #          arrowprops=dict(arrowstyle='<->'))
            #
            # ax.legend()
            # plt.show()
