import numpy as np
from astropy.io import fits
from pathlib import Path
from astro.papers.SHOC579_project.SHOC579_methods import VoxelPlotter, line_regions
import lime as lm
from matplotlib import pyplot as plt

obs_conf = lm.load_cfg(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\obsConf.ini')
files_data = obs_conf['data_location_windows']

data_folder = Path(files_data['data_folder'])
file_list = files_data['file_list']

percentil_array = obs_conf['sample_data']['percentil_array']
z_objs = obs_conf['sample_data']['z_array']
idx_j, idx_i = 37, 37

for i, file_name, in enumerate(file_list):

    file_address = data_folder/file_name
    with fits.open(file_address) as hdul:
        hdr = hdul['FLUX'].header
        wave = hdul['WAVE'].data
        flux = hdul['FLUX'].data
        err = 1/np.sqrt(hdul['IVAR'].data)
        pixMask = hdul['MASK'].data
    wave_rest = wave / (1+float(z_objs[i]))

    fig, ax = plt.subplots()
    ax.plot(wave, flux[:, idx_j, idx_i], label=f'Flux {idx_j}-{idx_i}')
    ax.plot(wave, err[:, idx_j, idx_i], label=f'Err {idx_j}-{idx_i}')
    ax.legend()
    ax.update({'xlabel': 'Wavelength', 'ylabel': 'Flux', 'title': 'Gaussian fitting'})
    plt.show()

    idcs_Halpha = np.searchsorted(wave, line_regions['H1_6563A'])
    Halpha_flux = flux[idcs_Halpha[0]:idcs_Halpha[1], :, :].sum(axis=0)
    levelContoursHalpha = np.nanpercentile(Halpha_flux, percentil_array)

    for lineLabel, lineLimits in line_regions.items():

        if lineLabel == 'O3_4363A':

            idcs_line = np.searchsorted(wave, lineLimits)
            lineMap = flux[idcs_line[0]:idcs_line[1], :, :].sum(axis=0)
            ion, wavelength, latexLabel = lm.label_decomposition(lineLabel, scalar_output=True)
            levelContours = np.nanpercentile(lineMap, percentil_array)

            plotConf = {'image': {'xlabel': r'RA', 'ylabel': r'DEC', 'title': f'MANGA SHOC579 {lineLabel}'}}
            VoxelPlotter(wave_rest, flux, Halpha_flux, (idx_j, idx_i), image_fg=lineMap,
                                   flux_levels=levelContours, ax_user_conf=plotConf, header=hdr)
