import numpy as np
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from src.specsiser.physical_model.line_tools import EmissionFitting, gauss_func

lm = sr.EmissionFitting()

STANDARD_PLOT = {'figure.figsize': (14, 14), 'axes.titlesize': 14, 'axes.labelsize': 14, 'legend.fontsize': 12,
                 'xtick.labelsize': 12, 'ytick.labelsize': 12}

STANDARD_AXES = {'xlabel': r'Wavelength $(\AA)$', 'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})$'}

rcParams.update(STANDARD_PLOT)

nrows, ncols = 5, 5

fig, ax = plt.subplots(nrows=nrows, ncols=ncols)
axesList = ax.flatten()

O3_4343A_params = np.array([1.030e-16, 4363.0, 2.35])

O3_4343A_cont = 1.05e-17
O3_4343A_label = r'$[OIII]4363\AA$'
mean_values = np.array([1.030e-16, 4363.0, 2.35])
ampTrue, muTrue, sigmaTrue = mean_values
stderr_values = np.array([ampTrue*0.10, muTrue*0.02, sigmaTrue*0.02])
wave_regions = np.array([4325, 4350, 4350, 4373, 4385, 4410])
areaTrue = np.sqrt(2 * np.pi * sigmaTrue ** 2) * ampTrue

lambda_res = np.array([0.50, 1, 1.5, 2.0, 2.5])
noise_res = np.array([0, 1.10, 1.25, 1.5, 2.0])

bbox_props = dict(boxstyle="round,pad=0.5", fc="w", ec="k", lw=1)

for i_colum in range(ncols):
    for i_row in range(nrows):

        wave = np.arange(4300, 4450, lambda_res[i_colum])
        flux_gauss = gauss_func((wave, O3_4343A_cont), ampTrue, muTrue, sigmaTrue)
        noise_level = np.random.normal(loc=O3_4343A_cont, scale=O3_4343A_cont*noise_res[i_row], size=wave.size)
        flux = flux_gauss + noise_level

        print(i_colum, i_row)

        # Measure the line flux
        lm = sr.EmissionFitting(wave, flux)

        idcsLines, idcsContinua = lm.define_masks(wave_regions)
        lm.line_properties(idcsLines, idcsContinua, bootstrap_size=500)
        lm.gauss_mcfit(idcsLines, idcsContinua, bootstrap_size=500)

        line_wave = wave[idcsLines]
        resampleWaveLine = np.linspace(line_wave[0], line_wave[-1], 100)
        resampleFluxCont = resampleWaveLine * lm.m_continuum + lm.n_continuum
        resampleGaussian = gauss_func((resampleWaveLine, resampleFluxCont), *lm.p1)
        idcsContinuumLeft = (wave_regions[0] <= lm.wave) & (lm.wave <= wave_regions[1])
        idcsContinuumRight = (wave_regions[4] <= lm.wave) & (lm.wave <= wave_regions[5])

        # Plot wavelength
        ax[i_row, i_colum].step(wave, flux)
        ax[i_row, i_colum].plot(resampleWaveLine, resampleGaussian)
        ax[i_row, i_colum].fill_between(wave[idcsLines], 0, flux[idcsLines], facecolor='tab:blue', step="pre", alpha=0.2)
        ax[i_row, i_colum].fill_between(wave[idcsContinuumLeft], 0, flux[idcsContinuumLeft], facecolor='tab:orange', step="pre", alpha=0.4)
        ax[i_row, i_colum].fill_between(wave[idcsContinuumRight], 0, flux[idcsContinuumRight], facecolor='tab:orange', step="pre", alpha=0.4)


        # Plot format
        line_ref = O3_4343A_label
        noise_ref = f'$\Delta F_{{\lambda}} = {lambda_res[i_row]}\cdot F_{{cont}}$'
        resolution_ref = f'$\Delta\lambda = {lambda_res[i_colum]}\AA$'
        ax[i_row, i_colum].yaxis.set_major_locator(plt.NullLocator())
        ax[i_row, i_colum].xaxis.set_major_locator(plt.NullLocator())
        ax[i_row, i_colum].update({'title': f'{line_ref}, {resolution_ref}, {noise_ref}'})
        flux_intg_comp = f'{lm.lineIntgFlux/areaTrue:0.2f}$\pm${lm.lineIntgErr/areaTrue:0.2f}'
        flux_gauss_comp = f'{lm.lineGaussFlux/areaTrue:0.2f}$\pm${lm.lineGaussErr/areaTrue:0.2f}'
        box_text = f'$F_{{intg}}$ = {flux_intg_comp}$\cdot F_{{true}}$\n' \
                   f'$F_{{gauss}}$ = {flux_gauss_comp}$\cdot F_{{true}}%$'
        ax[i_row, i_colum].annotate(box_text, xy=(0.6, 0.8), xycoords='axes fraction', xytext=None, bbox=bbox_props,
                                    verticalalignment='center')

plt.tight_layout()
plt.show()


