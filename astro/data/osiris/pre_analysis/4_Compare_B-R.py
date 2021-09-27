import numpy as np
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, gridspec
import astropy.units as u
# from src.specsiser.physical_model.line_tools import STANDARD_PLOT, STANDARD_AXES
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from src.specsiser.print.plot import STANDARD_PLOT

STANDARD_AXES = {'xlabel': r'Wavelength $(\AA)$', 'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})$'}

# Declare data and files location
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)
dataFolder = Path(obsData['file_information']['data_folder'])

objList = obsData['file_information']['object_list']
fileList = obsData['file_information']['files_list']
# objList_B = obsData['file_information']['objectB_list']
# fileList_B = obsData['file_information']['filesB_list']
# objList_R = obsData['file_information']['objectR_list']
# fileList_R = obsData['file_information']['filesR_list']
z_objs = obsData['sample_data']['z_array']
idx_band = int(obsData['file_information']['band_flux'])

for i, obj in enumerate(objList):

    fileList_B = dataFolder/fileList[i].replace('_BR', '_B')
    fileList_R = dataFolder/fileList[i].replace('_BR', '_R')

    file_BR, file_blue, file_red = dataFolder/fileList[i], dataFolder/fileList_B, dataFolder/fileList_R

    # Set and crop the wavelength
    rest_wave, flux, header = sr.import_fits_data(file_BR, instrument='OSIRIS')
    rest_waveB, fluxB, headerB = sr.import_fits_data(file_blue, instrument='OSIRIS')
    rest_waveR, fluxR, headerR = sr.import_fits_data(file_red, instrument='OSIRIS')

    z = z_objs[i]
    rest_wave, rest_waveB, rest_waveR = rest_wave/(1+z), rest_waveB/(1+z), rest_waveR/(1+z)

    # Plot Configuration
    defaultConf = STANDARD_PLOT.copy()
    rcParams.update(defaultConf)

    # Plot the spectra
    fig = plt.figure(figsize=(16,9))
    spec = gridspec.GridSpec(ncols=2, nrows=1, width_ratios=[16, 9])
    ax0 = fig.add_subplot(spec[0])
    ax0.step(rest_waveB, fluxB[idx_band][0], label='Blue arm ext=1', color='tab:blue')
    ax0.step(rest_waveR, fluxR[idx_band][0], label='Red arm ext=1', color='tab:red')
    ax0.step(rest_wave, flux, '--', label='Combined spectrum', color='tab:purple')

    idx_insetB = np.searchsorted(rest_waveB, (6200, 6800))
    idx_insetR = np.searchsorted(rest_waveR, (6200, 6800))
    idx_inset = np.searchsorted(rest_wave, (6200, 6800))

    ax1 = fig.add_subplot(spec[1])
    ax1.step(rest_waveB[idx_insetB[0]:idx_insetB[1]], fluxB[idx_band][0][idx_insetB[0]:idx_insetB[1]], label='Blue arm ext=1', color='tab:blue')
    ax1.step(rest_waveR[idx_insetR[0]:idx_insetR[1]], fluxR[idx_band][0][idx_insetR[0]:idx_insetR[1]], label='Red arm ext=1', color='tab:red')
    ax1.step(rest_wave[idx_inset[0]:idx_inset[1]], flux[idx_inset[0]:idx_inset[1]], '--', label='Combined spectrum', color='tab:purple')

    # axins = zoomed_inset_axes(ax, zoom=1.5, loc=8)
    # axins.step(rest_waveB[idx_insetB[0]:idx_insetB[1]], fluxB[idx_band][0][idx_insetB[0]:idx_insetB[1]], color='tab:blue')
    # axins.step(rest_waveR[idx_insetR[0]:idx_insetR[1]], fluxR[idx_band][0][idx_insetR[0]:idx_insetR[1]], color='tab:red')
    # axins.set_xlim(6200, 6800)
    # # axins.set_yscale('log')
    # # axins.set_ylim(-0.2e-16, 2.0e-16)
    # mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")

    # axins.get_xaxis().set_visible(False)
    # axins.get_yaxis().set_visible(False)

    ax0.set_yscale('log')
    ax0.update(STANDARD_AXES)
    ax0.set_title(f'Galaxy {obj}')
    ax0.legend()
    ax0.update(STANDARD_AXES)

    ax1.set_yscale('log')
    ax1.update(STANDARD_AXES)
    ax1.set_title(r'Galaxy {} $H\alpha$ region'.format(obj))
    ax1.legend()
    ax1.update(STANDARD_AXES)
    fig.tight_layout()

    plotAddress = dataFolder/fileList[i].replace('.fits', '_armFluxComparison.png')
    # plt.savefig(plotAddress, dpi=200, bbox_inches='tight')

    plt.show()

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