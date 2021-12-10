import src.specsiser as sr
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt, rcParams

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']
idx_band = int(obsData['file_information']['band_flux'])
z_objs = obsData['sample_data']['z_array']
flux_norm = obsData['sample_data']['norm_flux']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
plots_folder = Path(obsData['file_information']['images_folder'])

ext = 'BR'
cycle = 'it3'

papaderos_fittings_folder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/Papaderos_Full/6March2021/BCall_z03')
# papaderos_fittings_folder = Path('/home/vital/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/Papaderos_Full/6March2021/BCall_z03')



STANDARD_PLOT = {'figure.figsize': (10, 10),
                 'axes.titlesize': 14,
                 'axes.labelsize': 14,
                 'legend.fontsize': 14,
                 'xtick.labelsize': 12,
                 'ytick.labelsize': 12}
rcParams.update(STANDARD_PLOT)

# Reading spectrum
for i, obj in enumerate(objList):

    print(obj)
    if i == 5:

        # Declare input files
        objFolder = resultsFolder / f'{obj}'
        results_file = objFolder / f'{obj}_{ext}_measurements.txt'
        fits_file = dataFolder / f'{obj}_{ext}.fits'
        ext_fit = '1D'
        fado_file = f'{obj}_FD.cxt.FADO_{ext_fit}.fits'
        previousCycle = cycle.replace('3', '2')
        stellarFluxFile = objFolder / f'{obj}_{ext}_stellarFlux_{previousCycle}.txt'
        nebCompFile = objFolder/f'{obj}_{ext}_nebFlux_{previousCycle}.txt'

        # Load Fado data
        data, header = sr.import_fits_data(papaderos_fittings_folder/fado_file, instrument=None, frame_idx=0)
        wave_fado = np.arange(start=int(header['LAMBDA_I']), stop=int(header['LAMBDA_F'])+1)
        normF_fado = np.power(10, header["FLUXUNIT"]) * header["GALSNORM"]

        # Load observed spectrm
        wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
        wave_star, flux_star = np.loadtxt(stellarFluxFile, unpack=True)
        lm = sr.LineMesurer(wave, flux, redshift=z_objs[i], normFlux=normF_fado, crop_waves=(wmin_array[i], wmax_array[i]))
        wave_neb, flux_neb = np.loadtxt(nebCompFile, unpack=True)

        # Norm flux calculation
        normRange = (4125, 4300)
        idcs_normObs = (lm.wave > normRange[0]) & (lm.wave < normRange[1])
        idcs_normFado = (wave_fado > normRange[0]) & (wave_fado < normRange[1])
        normCoeff = np.mean(lm.flux[idcs_normObs])/np.mean(data[0][idcs_normFado])

        print(normCoeff)
        # fig, ax = plt.subplots(pdi=300)
        fig = plt.figure(dpi=300)
        ax = fig.add_subplot()

        ax.plot(lm.wave_rest, lm.flux/normCoeff*(1+lm.redshift), label='Observed spectrum')
        # ax.plot(wave_fado, data[0] * 1.0, label='Input spectrum Fado')
        # ax.plot(np.arange(len(data[2])), data[2], label='Masked pixels different than zero')

        ax.plot(wave_fado, data[3] * (1*1.03), label=f'FADO fit', linestyle='--')
        # ax.plot(wave_star, flux_star, label=f'STARLIGHT fit')
        ax.plot(wave_star, (flux_star + flux_neb)/normF_fado/normCoeff, label=f'STARLIGHT + nebular continuum fit', linestyle=':')
        ax.set_yscale('log')
        ax.set_xlim(3550, 4200)
        ax.set_ylim(1, 7)

        ax.legend()
        ax.update({'xlabel': r'Wavelength $(\AA)$', 'ylabel': r'Flux $(erg\,cm^{-2}s^{-1}\AA^{-1}) \cdot 1.132 \times 10^{-17}$'})
        plt.tight_layout()
        plt.savefig(plots_folder/'SSP_synthesis_comparison_300DPI.png')
        # plt.show()

