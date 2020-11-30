import numpy as np
from pathlib import Path
import src.specsiser as sr
import atpy
from src.specsiser.physical_model.starContinuum_functions import SSPsynthesizer, computeSSP_galaxy_mass
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astro.ext_lib.starlight.plotstarlightfits import plot_fits_and_SFH
import astro.ext_lib.starlight
from matplotlib import gridspec
import matplotlib
from matplotlib.ticker import ScalarFormatter


def sortHistBars(lpatch):
    x = np.array([p._x for p in lpatch])
    y = np.array([p._height for p in lpatch])
    zor = np.lexsort((y, x))
    d = {}
    for i, (ix, iy) in enumerate(zip(x[zor], y[zor])):
        d[(ix, iy)] = -i
    [p.set_zorder(d[(p._x, p._height)]) for p in lpatch]


objList = ['gp030321', 'gp101157', 'gp121903']
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList=objList, group_variables=False)
starlightFolder = Path(obsData['SSP_synthesis']['starlight_folder'])

fileList = obsData['file_information']['files_list']
dataFolder = Path(obsData['file_information']['data_folder'])
outputFolder = dataFolder/'flux_analysis'

objList_B = obsData['file_information']['objectB_list']
fileList_B = obsData['file_information']['filesB_list']
objList_R = obsData['file_information']['objectR_list']
fileList_R = obsData['file_information']['filesR_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])

counter = 0
for i, obj in enumerate(objList):

    z = z_objs[i]
    cHbeta = obsData[obj]['cHbeta']

    for ext in ('_BR', '_B'):

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        lineLog_file = outputFolder/f'{obj}{ext}_linesLog.txt'
        nebFluxNoNebCompFile = outputFolder/f'{obj}{ext}_obs_RemoveNebularComp.txt'
        nebCompFile = outputFolder/f'{obj}{ext}_NebFlux.txt'
        starlightOutput = starlightFolder/'Output'/f'{obj}{ext}.slOutput'
        objGasSpectrumFile = outputFolder/f'{obj}{ext}_gasSpectrum.txt'
        LightFracPlotFile = outputFolder/f'{obj}{ext}_SSP_LightFrac.png'
        stellarPlotFile = outputFolder/f'{obj}{ext}_stellarFit.png'
        specCompPlot = outputFolder/f'{obj}{ext}_specComponents.png'

        # Set wavelength and flux
        print(f'\n-- Treating {counter} :{obj}{ext}.fits')
        wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        wave_rest = wave / (1 + z)
        idx_wave = (wave_rest >= wmin_array[i]) & (wave_rest <= wmax_array[i])
        flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array

        obsWave, obsFlux = wave_rest[idx_wave], flux[idx_wave]
        nebWave, nebFlux = np.loadtxt(nebCompFile, unpack=True)
        obsNoNebWave, obsNoNebFlux = np.loadtxt(nebFluxNoNebCompFile, unpack=True)
        sw = SSPsynthesizer()
        stellarWave, inFlux, stellarFlux, fit_output = sw.load_starlight_output(starlightOutput)
        tsl = atpy.TableSet(str(starlightOutput), type='starlightv4')
        psfh = plot_fits_and_SFH(tsl)
        psfh.plot_fig_starlight()

        # collapse = False
        # cumsum = False
        # zcolor = [item['color'] for item in matplotlib.rcParams['axes.prop_cycle']]
        # fig = plt.figure(figsize=(10, 8))
        # gs = gridspec.GridSpec(2, 1)
        # gs.update(left=0.05, bottom=0.05, right=0.98, hspace=0.0)
        # ageBase = np.unique(tsl.population['popage_base'])
        # zBase = np.unique(tsl.population['popZ_base'])
        # lageBase = np.log10(ageBase)
        # log_age = np.log10(tsl.population['popage_base'])
        # # Bin width
        # d = log_age[1:] - log_age[:-1]
        # Bwidth = np.min(d[d > 0])
        # # Bwidth = np.append(np.diff(lageBase), np.diff(lageBase)[-1])
        # axl = plt.subplot(gs[0])  # Light
        # axm = plt.subplot(gs[1])  # Mass
        # patch_light = []
        # patch_mass = []
        # light_cum = np.zeros(lageBase.shape)
        # mass_cum = np.zeros(lageBase.shape)
        # for z, zcol in zip(zBase, zcolor):
        #     idz = tsl.population['popZ_base'] == z
        #     lbottom = light_cum if cumsum else None
        #     mbottom = mass_cum if cumsum else None
        #     if not collapse:
        #         pl = axl.bar(log_age[idz], tsl.population['popx'][idz], label='Z = %.4f' % z, width=Bwidth, color=zcol,
        #                      alpha=0.9, bottom=lbottom, align='center')
        #     patch_light.extend(pl)
        #     pm = axm.bar(log_age[idz], tsl.population['popmu_cor'][idz], label='Z = %.4f' % z, width=Bwidth, color=zcol,
        #                  alpha=0.9, bottom=mbottom, align='center')
        #     patch_mass.extend(pm)
        #     light_cum += tsl.population['popx'][idz]
        #     mass_cum += tsl.population['popmu_cor'][idz]
        # if collapse:
        #     pl = axl.bar(lageBase, light_cum, width=Bwidth, alpha=0.9)
        #     patch_light.extend(pl)
        #     pm = axm.bar(lageBase, mass_cum, width=Bwidth, alpha=0.9)
        #     patch_mass.extend(pm)
        # if not cumsum:
        #     sortHistBars(patch_light)
        #     sortHistBars(patch_mass)
        #
        # axl.set_ylabel('Light Fraction', fontsize=15)
        # axm.set_ylabel('Mass Fraction', fontsize=15)
        # axm.set_xlabel('log(Age)', fontsize=15)
        # if not collapse:
        #     axl.legend()
        #
        # axl.minorticks_on()
        # axl.tick_params(axis='x', which='major', labelbottom='off')
        # axm.minorticks_on()
        # axl.set_xlim(5.8, 10.5)
        # axm.set_xlim(5.8, 10.5)
        # axm.set_ylim(0., 90.)
        # axm.set_yscale('symlog', subsy=np.arange(1, 10))
        # axm.yaxis.set_major_formatter(ScalarFormatter())
        # fig.suptitle(tsl.keywords['arq_spec'], fontsize=20)
        plt.show()