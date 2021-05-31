import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, colors, cm, ticker, gridspec
from matplotlib.cbook import get_sample_data
from matplotlib.offsetbox   import OffsetImage, AnnotationBbox
from astro.data.muse.common_methods import STANDARD_AXES, DARK_PLOT, background_color, foreground_color

conf_file_address = Path('/home/vital/PycharmProjects/vital_tests/astro/papers/gtc_greenpeas/gtc_greenpeas_data.ini')
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])
plots_folder = Path('D:\Dropbox\Astrophysics\Papers\gtc_greenpeas\Greenpeas_osiris_manuscript_v0\images')

STANDARD_PLOT = {'figure.figsize': (7, 20),
                 'axes.titlesize': 14,
                 'axes.labelsize': 14,
                 'legend.fontsize': 12,
                 'xtick.labelsize': 12,
                 'ytick.labelsize': 12}
rcParams.update(STANDARD_PLOT)
format_image = {}
format_image['OffsetImage'] = {'zoom': 0.4}
format_image['AnnotationBbox'] = {'xy': (0.895, 0.65),
                                  'xybox': (0., 0.),
                                  'xycoords': 'axes fraction',
                                  'boxcoords': "offset points", "pad": 0.1}
counter = 0

for i, obj in enumerate(['gp030321', 'gp101157', 'gp121903']):

    z = z_objs[i]
    wmin, wmax = wmin_array[i], wmax_array[i]
    fit_conf = obsData[f'{obj}_line_fitting']

    for ext in ['_BR']:

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        objFolder = resultsFolder/f'{obj}'
        objMask = objFolder/f'{obj}{ext}_mask.txt'
        lineLog_file, lineGrid_file = objFolder/f'{obj}{ext}_linesLog.txt', objFolder/f'{obj}{ext}_lineGrid.png'
        pdfTableFile, txtTableFile = objFolder/f'{obj}{ext}_linesTable', objFolder/f'{obj}{ext}_linesTable.txt'
        image_obj = objFolder/f'{obj}_SDSS.png'
        print(f'\n-- Treating {i} :{obj}{ext}.fits')

        # wspace = 0, hspace = 0
        if obj in ['gp121903']:#['gp030321', 'gp101157', 'gp121903']:

            # Load spectrum
            wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
            flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
            lm = sr.LineMesurer(wave, flux, redshift=z_objs[i], normFlux=flux_norm, crop_waves=(wmin_array[i], wmax_array[i]))

            defaultConf = DARK_PLOT.copy()
            rcParams.update(defaultConf)

            fig = plt.figure(figsize=(16, 9))
            ax = fig.add_subplot()
            ax.step(lm.wave, lm.flux, color=foreground_color)
            # ax.set_yscale('log')
            ax.update({'xlabel': r'Wavelength $(\AA)$',
                       'ylabel': r'Flux $(10^{14}\,erg\,cm^{-2} s^{-1} \AA^{-1})$'})
            plt.tight_layout()
            plt.show()
