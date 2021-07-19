import numpy as np
from astropy.io import fits
from pathlib import Path
import src.specsiser as sr

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

verbose = False
52.963
17.601

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
    if verbose:
        lm.plot_spectrum(continuumFlux=lm.errFlux)

    # lm.plot_spectrum(obsLinesTable=obsLinesTable, matchedLinesDF=maskLinesDF, specLabel=f'{objName}')
    mask_local = data_folder/f'{objName}_mask'
    mask_local_df = sr.lineslogFile_to_DF(mask_local)
    if verbose:
        lm.plot_line_mask_selection(mask_local_df, mask_local, logscale=False)

    obsLines = mask_local_df.index.values
    for j, lineLabel in enumerate(obsLines):
        if lineLabel in ['O2_3726A_b', 'S2_6716A_b']:
            wave_regions = mask_local_df.loc[lineLabel, 'w1':'w6'].values
            lm.fit_from_wavelengths(lineLabel, wave_regions, user_conf=profile_conf)
            lm.print_results(show_plot=True, show_fit_report=True, log_scale=False)

    # Save the results
    lm.save_lineslog(lm.linesDF, lineslog_file)
    lm.table_fluxes(lm.linesDF, lineslog_table)
