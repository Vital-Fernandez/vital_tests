import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, gridspec
import astropy.units as u
from src.specsiser.components.line_tools import STANDARD_PLOT, STANDARD_AXES
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset

COLUMNS_TO_CLEAR = ['ion', 'intg_flux', 'intg_err', 'gauss_flux', 'gauss_err', 'eqw', 'eqw_err', 'm_cont',
                    'n_cont', 'std_cont', 'amp', 'mu', 'sigma', 'amp_err', 'mu_err', 'sigma_err',
                    'pynebCode', 'pynebLabel', 'lineType', 'blended_label', 'comments', 'cont']

# Declare data and files location
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)
objList = obsData['file_information']['object_list']
obsData = sr.loadConfData(conf_file_address, objList=objList, group_variables=False)

dataFolder = Path(obsData['file_information']['data_folder'])
objList = obsData['file_information']['object_list']
fileList = obsData['file_information']['files_list']
objList_B = obsData['file_information']['objectB_list']
fileList_B = obsData['file_information']['filesB_list']
objList_R = obsData['file_information']['objectR_list']
fileList_R = obsData['file_information']['filesR_list']
z_objs = obsData['sample_data']['z_array']
idx_band = int(obsData['file_information']['band_flux'])
wmin_array, wmax_array = obsData['sample_data']['wmin_array'], obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']

counter = 0
for i, obj in enumerate(objList):

    z = z_objs[i]
    # objMask = dataFolder/f'{obj}_BR_masks.txt'
    # maskDF = pd.read_csv(objMask, delim_whitespace=True, header=0, index_col=0)
    fit_conf = obsData[f'{obj}_line_fitting']

    for ext in ('_BR', '_B', '_R'):

            fits_address = dataFolder/f'{obj}{ext}.fits'
            lineLog_address = dataFolder/f'{obj}{ext}_linesLog_v0.txt'
            outputMask_file = dataFolder/'flux_analysis'/f'{obj}{ext}_mask.txt'
            linesLog_0_DF = pd.read_csv(lineLog_address, delim_whitespace=True, header=0, index_col=0)

            # Clear fit columns
            for column in COLUMNS_TO_CLEAR:
                if column in linesLog_0_DF.columns:
                    linesLog_0_DF.drop(column, axis=1, inplace=True)

            with open(outputMask_file, 'wb') as output_file:
                string_DF = linesLog_0_DF.to_string()
                output_file.write(string_DF.encode('UTF-8'))


            # Set and crop the wavelength
            print(f'-- Treating {counter} :{obj}{ext}.fits')
            wave, flux_array, header = sr.import_fits_data(fits_address, instrument='OSIRIS')
            wave_rest = wave/(1+z)

            if ext in ('_B', '_R'):
                flux = flux_array[idx_band][0]
            else:
                flux = flux_array

            # Define wave and flux ranges
            idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])

            # Load line measurer object
            lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave])

            # Find the lines
            lm.plot_spectrum(matchedLinesDF=linesLog_0_DF)

            # # Set and crop the wavelength
            # print(f'-- Treating {counter} :{obj}{ext}.fits')
            # wave, flux_array, header = sr.import_fits_data(fits_address, instrument='OSIRIS')
            # wave_rest = wave/(1+z)
            #
            # if ext in ('_B', '_R'):
            #     flux = flux_array[idx_band][0]
            # else:
            #     flux = flux_array

        #     # Define wave and flux ranges
        #     idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])
        #
        #     # Load line measurer object
        #     lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave])
        #
        #     # Find the lines
        #     norm_flux = lm.continuum_remover(noise_region)
        #     obsLinesTable = lm.line_finder(norm_flux, noiseWaveLim=noise_region, intLineThreshold=3)
        #     obsLinesDF = lm.match_lines(obsLinesTable, maskDF)
        #     lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=obsLinesDF)
        #
        #
        #     # # Get matched lines
        #     # idcsObsLines = (obsLinesDF.observation == 'detected')
        #     # obsLines = obsLinesDF.loc[idcsObsLines].index.values
        #     #
        #     # # Fit and check the regions
        #     # lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave])
        #     # for j, lineLabel in enumerate(obsLines):
        #     #     print(f'-- {lineLabel}:')
        #     #     wave_regions = obsLinesDF.loc[lineLabel, 'w1':'w6'].values
        #     #     lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=fit_conf)
        #     # lm.save_lineslog(lm.linesDF, lineLog_address)
        #
        #     # Plot checking model
        #     idcsObsLines = (obsLinesDF.observation == 'detected')
        #     lm.linesDF = obsLinesDF[idcsObsLines]
        #     lm.linesLogAddress = lineLog_address
        #     lm.plot_detected_lines(obsLinesDF[idcsObsLines], ncols=8)
        #
        # counter += 1

    # rest_wave, rest_waveB, rest_waveR = rest_wave/(1+z), rest_waveB/(1+z), rest_waveR/(1+z)
    #
    # # Plot Configuration
    # defaultConf = STANDARD_PLOT.copy()
    # rcParams.update(defaultConf)
    #
    # # Plot the spectra
    # fig = plt.figure(figsize=(16,9))
    # spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[16, 9])
    # ax0 = fig.add_subplot(spec[0])
    # ax0.step(rest_waveB, fluxB[idx_band][0], label='Blue arm ext=1', color='tab:blue')
    # ax0.step(rest_waveR, fluxR[idx_band][0], label='Red arm ext=1', color='tab:red')
    # ax0.step(rest_wave, flux, '--', label='Combined spectrum', color='tab:purple')
    #
    # idx_insetB = np.searchsorted(rest_waveB, (6200, 6800))
    # idx_insetR = np.searchsorted(rest_waveR, (6200, 6800))
    # idx_inset = np.searchsorted(rest_wave, (6200, 6800))
    #
    # ax1 = fig.add_subplot(spec[1])
    # ax1.step(rest_waveB[idx_insetB[0]:idx_insetB[1]], fluxB[idx_band][0][idx_insetB[0]:idx_insetB[1]], label='Blue arm ext=1', color='tab:blue')
    # ax1.step(rest_waveR[idx_insetR[0]:idx_insetR[1]], fluxR[idx_band][0][idx_insetR[0]:idx_insetR[1]], label='Red arm ext=1', color='tab:red')
    # ax1.step(rest_wave[idx_inset[0]:idx_inset[1]], flux[idx_inset[0]:idx_inset[1]], '--', label='Combined spectrum', color='tab:purple')
    #
    # # axins = zoomed_inset_axes(ax, zoom=1.5, loc=8)
    # # axins.step(rest_waveB[idx_insetB[0]:idx_insetB[1]], fluxB[idx_band][0][idx_insetB[0]:idx_insetB[1]], color='tab:blue')
    # # axins.step(rest_waveR[idx_insetR[0]:idx_insetR[1]], fluxR[idx_band][0][idx_insetR[0]:idx_insetR[1]], color='tab:red')
    # # axins.set_xlim(6200, 6800)
    # # # axins.set_yscale('log')
    # # # axins.set_ylim(-0.2e-16, 2.0e-16)
    # # mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")
    #
    # # axins.get_xaxis().set_visible(False)
    # # axins.get_yaxis().set_visible(False)
    #
    # ax0.set_yscale('log')
    # ax0.update(STANDARD_AXES)
    # ax0.set_title(f'Galaxy {obj}')
    # ax0.legend()
    # ax0.update(STANDARD_AXES)
    #
    # ax1.set_yscale('log')
    # ax1.update(STANDARD_AXES)
    # ax1.set_title(r'Galaxy {} $H\alpha$ region'.format(obj))
    # ax1.legend()
    # ax1.update(STANDARD_AXES)
    # fig.tight_layout()
    #
    # plotAddress = dataFolder/fileList[i].replace('.fits', '_armFluxComparison.png')
    # plt.savefig(plotAddress, dpi=200, bbox_inches='tight')

    # plt.show()

# fig, axes = plt.subplots(1,3, figsize = (12,4))
# x = np.arange(1,11)
# axes[0].plot(x, x**3, 'g',lw=2)
# axes[0].grid(True)
# axes[0].set_title('default grid')
# axes[1].plot(x, np.exp(x), 'r')
# axes[1].grid(color='b', ls = '-.', lw = 0.25)
# axes[1].set_title('custom grid')
# axes[2].plot(x,x)
# axes[2].set_title('no grid')
# fig.tight_layout()
# plt.show()