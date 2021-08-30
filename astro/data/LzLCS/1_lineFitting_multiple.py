import numpy as np
from astropy.io import fits
from pathlib import Path
import src.specsiser as sr

obsData = sr.loadConfData('./xshooter_LzLCS.ini')
data_folder = Path(obsData['data_location']['data_folder'])
results_folder = Path(obsData['data_location']['results_folder'])
# objfile_list = obsData['data_location']['objfile_list']
# sigmafile_list = obsData['data_location']['sigmafile_list']
# objRef_list = obsData['data_location']['ref_list']
maskfile = obsData['data_location']['generalMask']

obj_list = obsData['data_location']['obj_list']
ext_list = obsData['data_location']['ext_list']

wmin_array = obsData['sample_data']['w_min_array']
wmax_array = obsData['sample_data']['w_max_array']
norm_flux = obsData['sample_data']['norm_flux']
z_obj = obsData['sample_data']['z_obj_array']
profile_conf = obsData['line_fitting']

verbose = True

for i, objName in enumerate(obj_list):

    for j, ext in enumerate(ext_list):

        # Input data
        obj_folder = data_folder/f'sdss{objName}'
        spec_file = obj_folder/f'1dspectrum_{objName}_{ext}.fits'
        sigm_file = obj_folder/f'1dspectrum_{objName}_{ext}_sigma.fits'

        print(f'Treating {objName} {ext}, {z_obj[i]}')
        lineslog_file = results_folder/f'{objName}_{ext}_linesLog.txt'
        lineslog_table = results_folder/f'{objName}_{ext}_flux_table'
        obj_mask = data_folder/f'{objName}_{ext}_mask.txt'

        # Load inputs
        wave, flux, header = sr.import_fits_data(spec_file, instrument='xshooter', frame_idx=0)
        wave_sigma, sigma, header_sigma = sr.import_fits_data(sigm_file, instrument='xshooter', frame_idx=0)

        # Convert to angstroms
        wave = wave * 10 if objName != 'j131037' else wave
        wave_sigma = wave_sigma * 10 if objName != 'j131037' else wave_sigma

        trim_limits = [wmin_array[j], wmax_array[j]]
        lm = sr.LineMesurer(wave, flux, crop_waves=trim_limits, input_err=sigma, normFlux=norm_flux, redshift=z_obj[i])
        if verbose:
            lm.plot_spectrum(continuumFlux=lm.errFlux)

        # Find lines
        global_mask = data_folder/f'global_mask.txt'
        global_mask_df = sr.lineslogFile_to_DF(global_mask)
        norm_region = obsData['sample_data'][f'{ext}_norm_array']
        norm_spec = lm.continuum_remover(noiseRegionLims=norm_region)
        obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=norm_region, intLineThreshold=3)
        matchedDF = lm.match_lines(obsLinesTable, global_mask_df, find_line_borders=False)
        if verbose:
            lm.plot_spectrum(obsLinesTable=obsLinesTable, matchedLinesDF=matchedDF, specLabel=f'Emission line detection')

        # # lm.plot_spectrum(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF, specLabel=f'{objName}')
        # mask_local_df = sr.lineslogFile_to_DF(global_mask)
        # if verbose:
        #     lm.plot_line_mask_selection(mask_local_df, matchedDF, logscale=False)

        # obsLines = mask_local_df.index.values
        # for j, lineLabel in enumerate(obsLines):
        #     wave_regions = mask_local_df.loc[lineLabel, 'w1':'w6'].values
        #     lm.fit_from_wavelengths(lineLabel, wave_regions, user_conf=profile_conf)
        #     if verbose:
        #         lm.print_results(show_plot=True, show_fit_report=True, log_scale=False, frame='rest')

    # Save the results
    # lm.save_lineslog(lm.linesDF, lineslog_file)
    # lm.table_fluxes(lm.linesDF, lineslog_table)
