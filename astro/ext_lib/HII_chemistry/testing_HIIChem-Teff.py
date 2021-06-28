import numpy as np
import pandas as pd
import src.specsiser as sr
from physical_model.gasEmission_functions import gridInterpolatorFunction
from astro.papers.gtc_greenpeas.common_methods import epm_HII_CHI_mistry_Teff

conversion_dict = dict(O2_3726A_m='OII_3727',
                       Ne3_3869A='NeIII_3868',
                       O3_4363A='OIII_4363',
                       O3_4959A='OIII_4959',
                       O3_5007A='OIII_5007',
                       He1_4471A='HeI_4471',
                       He1_5876A='HeI_5876',
                       He2_4686A='HeII_4686',
                       N2_6584A='NII_6584',
                       S2_6716A='SII_6716',
                       S2_6731A='SII_6731',
                       S2_6716A_m='SII_6725')

# Set the true values
objParams = {'true_values': {'OH': 7.4,
                             'Teff': 60000.0,
                             'logU': -2.00,
                             'cHbeta': 0.0}}
errLines = 0.05
errAbund = 0.02

# errLines = 0.02
# errAbund = 0.01

# Load the ionization grids
gridLineDict, gridAxDict = sr.load_ionization_grid(log_scale=True)
lineInterpolator_dict = gridInterpolatorFunction(gridLineDict,
                                                 gridAxDict['logU'],
                                                 gridAxDict['Teff'],
                                                 gridAxDict['OH'],
                                                 interp_type='cube')

# Declare lines to simulate
lineLabels = np.array(['O2_3726A_m', 'He1_4471A', 'He2_4686A', 'O3_5007A', 'He1_5876A', 'S2_6716A_m'])
objParams['input_lines'] = lineLabels

# We use the default lines database to generate the synthetic emission lines log for this simulation
ion_array, wavelength_array, latexLabel_array = sr.label_decomposition(objParams['input_lines'])

# Define a pandas dataframe to contain the lines data
linesLogHeaders = ['wavelength', 'obsFlux', 'obsFluxErr', 'ion', 'blended_label', 'latexLabel']
objLinesDF = pd.DataFrame(index=lineLabels, columns=linesLogHeaders)
objLinesDF['ion'] = ion_array
objLinesDF['wavelength'] = wavelength_array
objLinesDF['latexLabel'] = latexLabel_array

# Declare extinction properties
objRed = sr.ExtinctionModel(Rv=3.1, red_curve='CCM89')

# Compute the reddening curve for the input emission lines
lineFlambdas = objRed.gasExtincParams(wave=objLinesDF.wavelength.values)

# Generate the line fluxes
lineLogFluxes = np.zeros(lineLabels.size)

coord_true = np.stack(([objParams['true_values']['logU']],
                       [objParams['true_values']['Teff']],
                       [objParams['true_values']['OH']]), axis=-1)

for i, lineLabel in enumerate(lineLabels):
    lineInt = lineInterpolator_dict[lineLabel](coord_true).eval()[0][0]
    print(i, lineLabel, np.power(10, lineInt), lineFlambdas[i])
    lineLogFluxes[i] = lineInt - objParams['true_values']['cHbeta'] * lineFlambdas[i]
lineFluxes = np.power(10, lineLogFluxes)

# Convert to a natural scale
objLinesDF['obsFlux'] = lineFluxes
objLinesDF['obsFluxErr'] = lineFluxes * errLines
objLinesDF.sort_values(by=['wavelength'], ascending=True, inplace=True)

# Epm style lines log
epm_log = pd.DataFrame()
epm_log.loc['0', '12logOH'] = objParams['true_values']['OH']
epm_log.loc['0', 'e12logOH'] = objParams['true_values']['OH'] * errAbund

# Recover the lines observed and label for HII-Chemistry
for sr_label, epm_label in conversion_dict.items():
    if sr_label in objLinesDF.index:
        epm_log.loc['0', epm_label] = objLinesDF.loc[sr_label, 'obsFlux']
        epm_log.loc['0', f'e{epm_label}'] = objLinesDF.loc[sr_label, 'obsFluxErr']

# Safe the dataframe
epmLineLog_file = 'D:/AstroModels/epm_lines_log_pp.txt'
with open(epmLineLog_file, 'wb') as output_file:
    string_DF = epm_log.to_string(index=False)
    output_file.write(string_DF.encode('UTF-8'))

# --------------------------- Run Epm fitting ----------------------------

# Sample data
HCm_folder = 'D:/Dropbox/Astrophysics/Tools/HCm-Teff_v5.01/'
grid_file = 'C17_bb_Teff_30-90_pp.dat'
epmLineLog_file = 'D:/AstroModels/epm_lines_log_pp.txt'
HII_Tef_fit_file = 'D:/AstroModels/HII_chem_output_pp.txt'

# HCm_folder = 'D:/Dropbox/Astrophysics/Tools/HCm-Teff_v5.01/'
# epmLineLog_file = 'D:/AstroModels/epm_lines_log_pp.txt'
# grid_file = 'C17_bb_Teff_30-90_sph.dat'
# HII_Tef_fit_file = 'D:/AstroModels/HII_chem_output_sph.txt'

# HCm_folder = 'D:/Dropbox/Astrophysics/Tools/HCm-Teff_v5.01/'
# epmLineLog_file = 'D:/AstroModels/epm_lines_log_pp_lowErr.txt'
# grid_file = 'C17_bb_Teff_30-90_pp.dat'
# HII_Tef_fit_file = 'D:/AstroModels/HII_chem_output_pp_lowErr.txt'

HCm_conf = dict(n=2000,
                sed=2,
                geo=1,
                inter=1,
                grid_file=f'{HCm_folder}/{grid_file}')

epm_HII_CHI_mistry_Teff(str(epmLineLog_file), HII_Tef_fit_file, **HCm_conf, HCm_folder=HCm_folder)

# Load fit data
output_headers = ['ID', 'O2Hb', 'eO2Hb', 'O3Hb', 'eO3Hb', '4471Hb', 'e4471Hb', '5876Hb', 'e5876Hb',
                  'He2Hb', 'eHe2Hb', 'S2Hb', 'eS2Hb', 'S3Hb', 'eS3Hb', 'O/H', 'eO/H', 'Teff', 'eTeff',
                  'logU', 'elogU']
HII_chem_array = np.loadtxt(HII_Tef_fit_file)
HII_chemDF = pd.DataFrame(columns=output_headers)
HII_chemDF.loc[0] = HII_chem_array

for item in HII_chemDF.columns.values:
    if (item is not 'ID') and (item[0] is not 'e'):
        errItem = f'e{item}'
        print(f'{item}: {HII_chemDF.loc[0, item]} +/- {HII_chemDF.loc[0, errItem]}')

# Compare observed fluxes with emission lines observed:
logU_fit = (HII_chemDF.loc[0, 'logU'], HII_chemDF.loc[0, 'elogU'])
Teff_fit = (HII_chemDF.loc[0, 'Teff'], HII_chemDF.loc[0, 'eTeff'])
OH_fit = (HII_chemDF.loc[0, 'O/H'], HII_chemDF.loc[0, 'eO/H'])

print(logU_fit)
print(Teff_fit)
print(OH_fit)

outputFluxes_dict = {}
MC_n = 1000

LogU_array = np.random.normal(logU_fit[0], logU_fit[1], size=MC_n)
Teff_vector = np.random.normal(Teff_fit[0], Teff_fit[1], size=MC_n)
OH_vector = np.random.normal(OH_fit[0], OH_fit[1], size=MC_n)
for i, lineLabel in enumerate(lineLabels):
    lineOutFluxes = np.zeros(MC_n)
    for i in np.arange(MC_n):
        coord_i = np.stack(([LogU_array[i]], [Teff_vector[i]], [OH_vector[i]]), axis=-1)
        lineOutFluxes[i] = lineInterpolator_dict[lineLabel](coord_i).eval()[0][0]
    outputFluxes_dict[lineLabel] = np.power(10, lineOutFluxes)

print('- Fitted fluxes')
for i, lineLabel in enumerate(lineLabels):
    print(lineLabel, np.mean(outputFluxes_dict[lineLabel]), np.std(outputFluxes_dict[lineLabel]))

# for line in
# for i in np.arange(MC_n):
#     lineOutFluxes = np.zeros(MC_n)
#     lineInt = lineInterpolator_dict[lineLabel](coord_true).eval()[0][0]
#     print(i, lineLabel, np.power(10, lineInt), lineFlambdas[i])
#     lineLogFluxes[i] = lineInt - objParams['true_values']['cHbeta'] * lineFlambdas[i]
# lineFluxes = np.power(10, lineLogFluxes)
