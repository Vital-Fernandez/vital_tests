import numpy as np
import pandas as pd
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from pathlib import Path
from lmfit.models import GaussianModel, LinearModel

# Import the observation data
obsData = sr.loadConfData('gtc_greenpeas_data.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['file_information']['data_folder'])
file_list = obsData['file_information']['files_list']
addressList = list(data_folder / file for file in file_list)

# Get the file
idx_file = 0
fits_address = addressList[idx_file]
linesLog_address = str(fits_address).replace('.fits', '_linesLog.txt')
linesDF = pd.read_csv(linesLog_address, delim_whitespace=True, header=0, index_col=0)
objRef = obsData['file_information']['object_list'][idx_file]
print(f'- Treating file: {obsData["file_information"]["files_list"][idx_file]}')

# Get fits data
norm_flux = obsData['sample_data']['norm_flux']
wmin_array, wmax_array = obsData['sample_data']['wmin_array'], obsData['sample_data']['wmax_array']
wave, flux, header = sr.import_fits_data(fits_address, instrument='OSIRIS')
wave_rest = wave / (1 + obsData['sample_data']['z_array'][idx_file])
idx_wave = (wave_rest >= wmin_array[idx_file]) & (wave_rest <= wmax_array[idx_file])

# Analyse the spectrum
lm = sr.LineMeasurer(wave_rest[idx_wave], flux[idx_wave] / norm_flux)
# lm.plot_spectrum_components()

# Measure line fluxes
lineLabel = 'H1_6563A'
wave_regions = linesDF.loc[lineLabel, 'w1':'w6'].values
lineBlendComp = linesDb.loc[lineLabel, 'blended']

idcsLinePeak, idcsContinua = lm.define_masks(wave_regions)
lm.line_properties(idcsLinePeak, idcsContinua, bootstrap_size=500)
lm.line_fitting(idcsLinePeak, idcsContinua, bootstrap_size=500)

idx_line = np.searchsorted(lm.wave, [wave_regions[0], wave_regions[-1]])
wave_region, flux_region = lm.wave[idx_line[0]:idx_line[1]], lm.flux[idx_line[0]:idx_line[1]]
lineComponents = obsData[objRef]['H1_6563A_b_components']

obsData[f'{objRef}_blended_lines']['cont_slope']['value'] = lm.m_continuum
obsData[f'{objRef}_blended_lines']['cont_intercept']['value'] = lm.n_continuum
obsData[f'{objRef}_blended_lines']['H1_6563A_amplitude']['value'] = lm.peakInt
obsData[f'{objRef}_blended_lines']['H1_6563A_center']['value'] = lm.peakWave

fitOutput = lm.composite_fitting(lineComponents, idcsLinePeak, idcsContinua, continuum_check=True,
                                 user_conf=obsData[f'{objRef}_blended_lines'], wide_comp='H1_6563Aw')
print(fitOutput.fit_report())

param_value, param_stderr = fitOutput.params['N2_6584A_center'].value, fitOutput.params['N2_6584A_center'].stderr

lm.plot_fit_components(fitOutput)

