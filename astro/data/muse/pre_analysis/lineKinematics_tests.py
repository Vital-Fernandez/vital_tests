import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import src.specsiser as sr
from pathlib import Path
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, red_corr_HalphaHbeta_ratio, default_linelog_types
from astropy.io import fits
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
import time

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

dict_errs = {}
dict_nan_values = {}

ref_flux_line = 'S3_6312A'
sulfur_bdry = int(obsData['Region_masks']['S3_direct_limit'])
hydrogen_bdry = int(obsData['Region_masks']['H1_direct_limit'])
verbose = True
outSide_file = Path('/home/vital/Astro-data/muse_voxel_170-170.txt')
outSide_df = Path('/home/vital/Astro-data/muse_voxel_170-170_df.txt')
outSide_kinTable = Path('/home/vital/Astro-data/muse_voxel_170-170_kinematics')
outSide_kin_df = Path('/home/vital/Astro-data/muse_voxel_170-170_KinDF.txt')
obj = objList[0]

# Data location
objFolder = resultsFolder/obj

mask_address = dataFolder/obj/f'{obj}_mask.txt'
mask_df = pd.read_csv(mask_address, delim_whitespace=True, header=0, index_col=0)

wave, flux, err_flux = np.loadtxt(outSide_file, unpack=True)

# cube_address = fitsFolder / fileList[0]
# wave, cube, header = sr.import_fits_data(cube_address, instrument='MUSE')
# flux_voxel = cube[:, 170, 170].data
# flux_err = cube[:, 170, 170].var


# # --------------------------------- Single line fitting -------------------------------------
# lm = sr.LineMesurer(wave, flux, input_err=err_flux, normFlux=norm_flux)
#
# lineLabel = 'H1_6563A_b'
# # wave_regions = mask_df.loc['H1_6563A', 'w1':'w6'].values
# wave_regions = np.array([6499.49, 6526.84, 6530., 6592., 6593.2, 6616.02])
#
# # line_conf = {'H1_6563A_b': 'H1_6563A-H1_6563A_w1-N2_6584A-N2_6548A'}
# line_conf = {'H1_6563A_b': 'H1_6563A-H1_6563A_w1-H1_6563A_w2-N2_6584A-N2_6548A'}
# # line_conf = {'H1_6563A_b': 'H1_6563A-H1_6563A_w1-N2_6584A-N2_6548A'}
#
# line_conf['N2_6548A_sigma'] = {'expr': 'N2_6584A_sigma'}
# line_conf['N2_6548A_amplitude'] = {'expr': 'N2_6584A_amplitude/2.94'}
# line_conf['H1_6563A_w1_sigma'] = {'expr': '>2*H1_6563A_sigma'}
# line_conf['H1_6563A_w1_sigma'] = {'expr': '>2*H1_6563A_sigma'}
# line_conf['H1_6563A_w2_sigma'] = {'expr': '>2*H1_6563A_w1_sigma'}
#
# lm.fit_from_wavelengths('H1_6563A_b', wave_regions, fit_conf=line_conf)
# # lm.print_results(show_fit_report=True)
# lm.plot_fit_components(lm.fit_output, logScale=True)
# lm.save_lineslog(lm.linesDF, file_address=outSide_df)

#--------------------------------- Multiple lines -------------------------------------
# lm = sr.LineMesurer(wave, flux, input_err=err_flux, normFlux=norm_flux)
# # lm.plot_spectrum_components(specLabel=f'{obj} voxel', log_scale=True)
# # Security check for pixels with nan values:
# idcs_nan = np.isnan(lm.flux)
#
# if idcs_nan.any():
#     Interpolation = interp1d(lm.wave[~idcs_nan], lm.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
#     lm.flux = Interpolation(lm.wave)
#     norm_spec = lm.continuum_remover(noise_region)
# else:
#     norm_spec = lm.continuum_remover(noise_region)
#
# # Identify the emission lines
# obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=noise_region, intLineThreshold=3)
# maskLinesDF = lm.match_lines(obsLinesTable, mask_df)
# # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF,
# #                             specLabel=f'{obj} voxel')
#
# line_conf = {'H1_6563A_b': 'H1_6563A-H1_6563A_w1-H1_6563A_w2-N2_6584A-N2_6548A'}
# line_conf['N2_6548A_sigma'] = {'expr': 'N2_6584A_sigma'}
# line_conf['N2_6548A_amplitude'] = {'expr': 'N2_6584A_amplitude/2.94'}
# line_conf['H1_6563A_w1_sigma'] = {'expr': '>2*H1_6563A_sigma'}
# line_conf['H1_6563A_w1_sigma'] = {'expr': '>2*H1_6563A_sigma'}
# line_conf['H1_6563A_w2_sigma'] = {'expr': '>2*H1_6563A_w1_sigma'}
# maskLinesDF.rename(index={'H1_6563A': 'H1_6563A_b'}, inplace=True)
# maskLinesDF.drop(index=['N2_6584A', 'N2_6548A'], inplace=True)
# maskLinesDF.loc['H1_6563A_b', 'w1':'w6'] = np.array([6499.49, 6526.84, 6530., 6592., 6593.2, 6616.02])
#
#
# idcsObsLines = (maskLinesDF.observation == 'detected')
#
# # Reset and measure the lines
# lm = sr.LineMesurer(wave, flux, input_err=err_flux, redshift=0, normFlux=norm_flux)
# obsLines = maskLinesDF.loc[idcsObsLines].index.values
# for j, lineLabel in enumerate(obsLines):
#     wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
#     lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=line_conf)
#     # lm.plot_fit_components(lm.fit_output, logScale=True)
# lm.save_lineslog(lm.linesDF, file_address=outSide_df)


# #--------------------------------- Make table -------------------------------------
# lm = sr.LineMesurer(wave, flux, input_err=err_flux, normFlux=norm_flux)
# linesDF = sr.lineslogFile_to_DF(outSide_df)
# lm.table_kinematics(linesDF, outSide_kinTable)

#--------------------------------- Copying kinematics -------------------------------------
lm = sr.LineMesurer(wave, flux, input_err=err_flux, normFlux=norm_flux)
# linesDF = sr.lineslogFile_to_DF(outSide_df)

obsLines = mask_df.index.to_list()
O3_lines = ['O3_5007A', 'O3_4959A', 'O3_4363A']
a = sorted(obsLines, key=lambda e: (O3_lines.index(e), e) if e in O3_lines else (len(O3_lines), e))
mask_df = mask_df.reindex(index=a)
mask_df.rename(index={'H1_6563A': 'H1_6563A_b'}, inplace=True)
mask_df.drop(index=['N2_6584A', 'N2_6548A'], inplace=True)
mask_df.loc['H1_6563A_b', 'w1':'w6'] = np.array([6499.49, 6526.84, 6530., 6592., 6593.2, 6616.02])

line_conf = {'H1_6563A_b': 'H1_6563A-H1_6563A_w1-H1_6563A_w2-N2_6584A-N2_6548A'}
line_conf['N2_6548A_sigma'] = {'expr': 'N2_6584A_sigma'}
line_conf['N2_6548A_kinem'] = 'N2_6584A'
line_conf['N2_6548A_amplitude'] = {'expr': 'N2_6584A_amplitude/2.94'}
line_conf['H1_6563A_w1_sigma'] = {'expr': '>2*H1_6563A_sigma'}
line_conf['H1_6563A_w1_sigma'] = {'expr': '>2*H1_6563A_sigma'}
line_conf['H1_6563A_w2_sigma'] = {'expr': '>2*H1_6563A_w1_sigma'}
line_conf['O3_4959A_kinem'] = 'O3_5007A'

for j, lineLabel in enumerate(mask_df.index.values):
    if lineLabel == 'H1_6563A_b':
        print(lineLabel, )
        wave_regions = mask_df.loc[lineLabel, 'w1':'w6'].values
        lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=line_conf)
        lm.print_results(show_fit_report=True)
        lm.plot_fit_components(lm.fit_output, logScale=True)
lm.save_lineslog(lm.linesDF, file_address=outSide_kin_df)
lm.table_kinematics(lm.linesDF, table_address=Path('/home/vital/Astro-data/muse_voxel_170-170_kinematics_IMPORT'))
#--------------------------------- Code failing warning -------------------------------------
# lm = sr.LineMesurer(wave, flux, input_err=err_flux, normFlux=norm_flux)
# # lm.plot_spectrum_components(specLabel=f'{obj} voxel', log_scale=True)
# # Security check for pixels with nan values:
# idcs_nan = np.isnan(lm.flux)
#
# if idcs_nan.any():
#     Interpolation = interp1d(lm.wave[~idcs_nan], lm.flux[~idcs_nan], kind='slinear', fill_value="extrapolate")
#     lm.flux = Interpolation(lm.wave)
#     norm_spec = lm.continuum_remover(noise_region)
# else:
#     norm_spec = lm.continuum_remover(noise_region)
#
# # Identify the emission lines
# obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=noise_region, intLineThreshold=3)
# maskLinesDF = lm.match_lines(obsLinesTable, mask_df)
# # lm.plot_spectrum_components(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF,
# #                             specLabel=f'{obj} voxel')
#
# line_conf = {'H1_6563A_b': 'H1_6563A-H1_6563A_w1-H1_6563A_w2-N2_6584A-N2_6548A'}
# line_conf['H1_6563A_w1_sigma'] = {'expr': '>2*H1_6563A_sigma'}
# line_conf['H1_6563A_w1_sigma'] = {'expr': '>2*H1_6563A_sigma'}
# line_conf['H1_6563A_w2_sigma'] = {'expr': '>2*H1_6563A_w1_sigma'}
# line_conf['N2_6548A_sigma'] = {'expr': 'N2_6584A_sigma'}
# line_conf['N2_6548A_amplitude'] = {'expr': 'N2_6584A_amplitude/2.94'}
# maskLinesDF.rename(index={'H1_6563A': 'H1_6563A_b'}, inplace=True)
# maskLinesDF.drop(index=['N2_6584A', 'N2_6548A'], inplace=True)
#
# idcsObsLines = (maskLinesDF.index == 'H1_6563A_b')
#
# # Reset and measure the lines
# lm = sr.LineMesurer(wave, flux, input_err=err_flux, redshift=0, normFlux=norm_flux)
# obsLines = maskLinesDF.loc[idcsObsLines].index.values
# for j, lineLabel in enumerate(obsLines):
#     wave_regions = maskLinesDF.loc[lineLabel, 'w1':'w6'].values
#     lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=line_conf)
#     # lm.plot_fit_components(lm.fit_output, logScale=True)
# lm.save_lineslog(lm.linesDF, file_address=outSide_df)