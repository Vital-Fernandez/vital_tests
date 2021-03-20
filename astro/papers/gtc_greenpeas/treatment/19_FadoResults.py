import src.specsiser as sr
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

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

ext = 'BR'
cycle = 'it3'

# papaderos_fittings_folder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/Papaderos_Full/out_FD')
# ssp_lib_list = ['ChabP94Z5N295', 'midChab94', 'Z3SalpP00', 'Z4P00']
#
# """
# _1D: One dimensional spectra contains:
# _EL: Emission-line results
# _ST: Statistics: By-products from FADO
# _DE: Differential evolution parameters
# """
#
# ext_fado_fits_list = ['1D', 'EL', 'ST', 'DE']
# params_per_ext = {'1D':['ARQ_BASE', 'ARQ_CONF', 'CHI2_VAL', 'CHI2_RED', 'NUMPARAM', 'NUM_BASE'],
#                   'EL':['NE_LINES', 'N_PARAMS', 'LAMBDA_0', 'GALSNORM', 'TELECTRO', 'DELECTRO', 'GEXTINCT', 'GEXTBDEV', 'GNEBULAR', 'GNEBBDEV'],
#                   'ST':['GEXTINCT', 'GNEBULAR', 'AVE_LAGE', 'DEV_LAGE', 'AVE_MAGE', 'DEV_MAGE', 'AVELLAGE', 'DEVLLAGE', 'AVELMAGE', 'DEVLMAGE',
#                         'AVE_LMET', 'DEV_LMET', 'AVE_MMET', 'DEV_MMET',
#                         'LOGPEAVE', 'LOGPEDEV', 'LOGPCAVE', 'LOGPCDEV',
#                         'LOGMCAVE', 'LOGMCDEV', 'LOGPEAVE', 'LOGPEDEV',
#                         'LOGMEAVE', 'LOGMEDEV'],
#                   'DE':[]}
#
# # Reading parameters from four libraries
# for i, obj in enumerate(objList):
#
#     # Declare input files
#     objFolder = resultsFolder / f'{obj}'
#     results_file = objFolder / f'{obj}_{ext}_measurements.txt'
#
#     # # Load the data
#     # results_dict = sr.loadConfData(results_file, group_variables=False)
#
#     for j, ssp_lib in enumerate(ssp_lib_list):
#
#         # Results container
#         fado_measurements = {}
#         for k, ext_fit in enumerate(ext_fado_fits_list):
#
#             fits_address = papaderos_fittings_folder/ssp_lib/f'{obj}_FD.cxt.FADO_{ext_fit}.fits'
#             data, header = sr.import_fits_data(fits_address, instrument=None, frame_idx=0)
#
#             for param in params_per_ext[ext_fit]:
#                 print(f'In {ssp_lib} fittig: {param} = {header[param]} ')
#                 fado_measurements[param] = header[param]
#
#         # Store the results
#         print(fado_measurements)
#         section_label = f'{ssp_lib}_FADO'
#         sr.parseConfDict(results_file, fado_measurements, section_label, clear_section=True)


papaderos_fittings_folder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/Papaderos_Full/6March2021/BCall_z03')

params_per_ext = {'1D': ['ARQ_BASE', 'ARQ_CONF', 'CHI2_VAL', 'CHI2_RED', 'NUMPARAM', 'NUM_BASE'],
                  'EL': ['NE_LINES', 'N_PARAMS', 'LAMBDA_0', 'GALSNORM', 'TELECTRO', 'DELECTRO'],
                  'ST': ['BST_LAGE', 'DEV_LAGE',
                         'BSTLLAGE', 'DEVLLAGE',
                         'BST_MAGE', 'DEV_MAGE',
                         'BSTLMAGE', 'DEVLMAGE',
                         # 'BSTLMMET', 'DEVLMMET',
                         'BST_LMET', 'DEV_LMET',
                         'BST_MMET', 'DEV_MMET',
                         'LOGMEBST', 'LOGMEDEV',
                         'LOGMCBST', 'LOGMCDEV',
                         'LOGPEBST', 'LOGPEDEV',
                         'LOGPCBST', 'LOGPCDEV',
                         'GEXTINCT', 'GEXTBDEV',
                         'GNEBULAR', 'GNEBBDEV']}

# Reading parameters from four libraries
for i, obj in enumerate(objList):

    # Declare input files
    objFolder = resultsFolder / f'{obj}'
    results_file = objFolder / f'{obj}_{ext}_measurements.txt'

    # Results container
    fado_measurements = {}

    for ext_fit, list_params in params_per_ext.items():

        ssp_lib = f'{obj}_FD.cxt.FADO_{ext_fit}.fits'
        fits_address = papaderos_fittings_folder/ssp_lib
        data, header = sr.import_fits_data(fits_address, instrument=None, frame_idx=0)

        for param in list_params:
            # print(f'{param} = {header[param]}')
            fado_measurements[param] = header[param]

        # for param in header:
        #     print(f'{param} = {header[param]}')

    # Store the results
    section_label = f'Large_FADO_fit'
    sr.parseConfDict(results_file, fado_measurements, section_label, clear_section=True)


# papaderos_fittings_folder = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/Papaderos_Full/6March2021/BCall_z03')
# ssp_lib_list = ['ChabP94Z5N295', 'midChab94', 'Z3SalpP00', 'Z4P00']
#
# # Reading spectrum
# for i, obj in enumerate(objList):
#
#     # Declare input files
#     objFolder = resultsFolder / f'{obj}'
#     results_file = objFolder / f'{obj}_{ext}_measurements.txt'
#     fits_file = dataFolder / f'{obj}_{ext}.fits'
#     ext_fit = '1D'
#     fado_file = f'{obj}_FD.cxt.FADO_{ext_fit}.fits'
#     previousCycle = cycle.replace('3', '2')
#     stellarFluxFile = objFolder / f'{obj}_{ext}_stellarFlux_{previousCycle}.txt'
#     nebCompFile = objFolder/f'{obj}_{ext}_nebFlux_{previousCycle}.txt'
#
#     # Load Fado data
#     data, header = sr.import_fits_data(papaderos_fittings_folder/fado_file, instrument=None, frame_idx=0)
#     wave_fado = np.arange(start=int(header['LAMBDA_I']), stop=int(header['LAMBDA_F'])+1)
#     normF_fado = np.power(10, header["FLUXUNIT"]) * header["GALSNORM"]
#
#     # Load observed spectrm
#     wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
#     flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
#     wave_star, flux_star = np.loadtxt(stellarFluxFile, unpack=True)
#     lm = sr.LineMesurer(wave, flux, redshift=z_objs[i], normFlux=normF_fado, crop_waves=(wmin_array[i], wmax_array[i]))
#     wave_neb, flux_neb = np.loadtxt(nebCompFile, unpack=True)
#
#     fig, ax = plt.subplots()
#     ax.plot(lm.wave, lm.flux, label='Observed spectrum')
#     ax.plot(wave_fado, data[0] * 1.0, label='Input spectrum Fado')
#     # ax.plot(np.arange(len(data[2])), data[2], label='Masked pixels different than zero')
#
#     ax.plot(wave_fado, data[3] * 1.0, label=f'FADO fit')
#     # ax.plot(wave_star, flux_star, label=f'STARLIGHT fit')
#     ax.plot(wave_star, (flux_star + flux_neb)/normF_fado, label=f'Stellar + nebular fit')
#     ax.set_yscale('log')
#     ax.legend()
#     ax.update({'xlabel': 'Wavelength', 'ylabel': 'Flux', 'title': 'FADO outputs'})
#     plt.show()
