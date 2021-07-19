import numpy as np
from astropy.io import fits
from pathlib import Path
from src.specsiser.data_printing import DARK_PLOT
import src.specsiser as sr
import pyneb as pn
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

verbose=True

Ne3 = pn.Atom('Ne', 3)
H1 = pn.RecAtom('H', 1)

temp, den = 10000, 100

Ne3_alone = Ne3.getEmissivity(temp, den, wave=3869)/H1.getEmissivity(temp, den, wave=4861)
Ne3_H7 = Ne3.getEmissivity(temp, den, wave=3970)/H1.getEmissivity(temp, den, wave=4861)

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
    lm.plot_spectrum(continuumFlux=lm.errFlux, plotConf=DARK_PLOT)