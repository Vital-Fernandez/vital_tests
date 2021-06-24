import numpy as np
import pandas as pd
from matplotlib import pyplot as plt, rcParams
from pathlib import Path
import src.specsiser as ss


linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path('D:/Dropbox/Astrophysics/Data/WHT-Ricardo/')
fileList = ['COMBINED_blue.0001.fits', 'combined_red.0001.fits']
addressList = list(data_folder/file for file in fileList)

# Get the fits file data
wave_blue, flux_blue, header_blue = ss.import_fits_data(addressList[0], 0)
wave_red, flux_red, header_red = ss.import_fits_data(addressList[1], 0)

# Redshift correction
redshift, blue_limits = 1.1735, (3600, 5100)
wave_blue = wave_blue/redshift
wave_red = wave_red/redshift

# Trim the observation
idcs_limit = (blue_limits[0] <= wave_blue) & (wave_blue <= blue_limits[1])
wave_blue, flux_blue = wave_blue[idcs_limit], flux_blue[idcs_limit]

lm = ss.EmissionFitting(wave_blue, flux_blue)

# Remove the continuum
flux_noContinuum = lm.continuum_remover(noiseWaveLim=(4150, 4300))

# Find lines
linesTable = lm.line_finder(flux_noContinuum, intLineThreshold=4, noiseWaveLim=(4150, 4300), verbose=False)

# Match lines
linesDb = lm.match_lines(linesTable, linesDb)

# Plot the complete spectrum
lm.plot_spectrum_components(lm.flux - flux_noContinuum, linesTable, linesDb)

# Plot the matched lines:
lm.plot_line_mask_selection(linesDb)

# Measure line fluxes
idcsObsLines = (linesDb.observation == 'detected')
obsLines = linesDb.loc[idcsObsLines].index.values

for i in np.arange(obsLines.size):

    lineLabel = obsLines[i]
    print(f'- {lineLabel}:')

    # Declare regions data
    wave_regions = linesDb.loc[lineLabel, 'w1':'w6'].values
    idcsLinePeak, idcsContinua = lm.define_masks(wave_regions)

    # Identify line regions
    lm.line_properties(idcsLinePeak, idcsContinua, bootstrap_size=500)

    # Perform gaussian fitting
    lm.gauss_mcfit(idcsLinePeak, idcsContinua, bootstrap_size=500)

    # Store results in database
    lm.results_to_database(lineLabel, linesDb)

# Save dataframe to text file
linesLogAddress = data_folder / fileList[0].replace('.fits', '_linesLog.txt')
lm.save_lineslog(linesDb.loc[idcsObsLines], linesLogAddress)

# Plot the matched lines:
lm.plot_line_mask_selection(linesDb)


# for label in ['H1_4861A', 'Ne3_3869A', 'H1_4341A', 'O3_5007A']:
#
#     waveLimits = linesDb.loc[label, 'w1':'w6'].values
#
#     idcsLine, idcsContinua = lm.define_masks(waveLimits)
#
#     # Identify line regions
#     lm.line_properties(idcsLine, idcsContinua, bootstrap_size=1000)
#
#     # Fit gaussian profile
#     lm.line_fitting(idcsLine, idcsContinua, bootstrap_size=1000)
#
#     # Pre-print
#     outputLine = f'Integrated: {lm.intg_flux}+/-{lm.intg_err} -- Gaussian: {lm.gauss_flux}+/-{lm.gauss_err}'
#     print(outputLine)

# fig, ax = plt.subplots()
# ax.plot(lm.wave[idcsContinua], lm.flux[idcsContinua], label='continuum')
# ax.plot(lm.wave[idcsLine], lm.flux[idcsLine], label='Line')
# ax.legend()
# ax.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'Gaussian fitting'})
# plt.show()



# idcsLineType = linesTable['line_type'] == 'emission'
# idcsLinePeak = np.array(linesTable[idcsLineType]['line_center_index'])
#
#
# fig, ax = plt.subplots()
# ax.plot(lm.wave, lm.flux, label='Blue arm')
# ax.plot(lm.wave, lm.flux-flux_noContinuum, label='Blue arm')
# ax.scatter(lm.wave[idcsLinePeak], lm.flux[idcsLinePeak], label='linePeaks', facecolors='tab:green')
# # ax.plot(wave_blue, flux_blue, label='Blue arm')
# # ax.plot(wave_red, flux_red, label='Red arm')
# ax.legend()
# ax.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'Gaussian fitting'})
# plt.show()
#
# print('hi')


