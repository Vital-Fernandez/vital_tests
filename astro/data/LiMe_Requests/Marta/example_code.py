import numpy as np
import lime




def HII_flux(wave, flux):

    # State the data files
    lineMaskFile = './SDSS_BR_mask.txt'

    # Selection of reference for the plots
    plots_frame = 'obs'

    # # Load configuration
    # cfgFile = './config_file.cfg'
    # obs_cfg = lime.load_cfg(cfgFile)
    # fit_cfg = obs_cfg['SDSS_line_fitting']

    # Load mask
    maskDF = lime.load_lines_log(lineMaskFile)

    z_obj = 0
    norm_flux_obj = 1e-17

    # Line name and its location mask in the rest frame
    lineLabel = 'H1_6563A_b'
    lineWaves = np.array([6438.03, 6508.66, 6535.10, 6600.95, 6627.70, 6661.82])

    # Define a spectrum object
    sdss_spec = lime.Spectrum(wave, flux, redshift=z_obj, norm_flux=norm_flux_obj)
    # sdss_spec.plot_spectrum()

    # Find lines
    peaks_table, matched_masks_DF = sdss_spec.match_line_mask(maskDF, [5600, 5850])
    # sdss_spec.plot_spectrum(peaks_table=peaks_table, matched_DF=matched_masks_DF, spec_label=f'SDSS spectrum')

    # Correct line region
    corrected_mask_file = '.SDSS_BR_mask_corrected.txt'
    lime.save_line_log(matched_masks_DF, corrected_mask_file)

    # Object line fitting configuration

    # Measure the emission lines
    for i, lineLabel in enumerate(matched_masks_DF.index.values):
        wave_regions = matched_masks_DF.loc[lineLabel, 'w1':'w6'].values
        sdss_spec.fit_from_wavelengths(lineLabel, wave_regions)  # , user_cfg=fit_cfg)
        # sdss_spec.display_results(fit_report=True, plot=True, log_scale=True, frame='obs')

    # Save the results
    sdss_spec.plot_line_grid(sdss_spec.log)
    lime.save_line_log(sdss_spec.log, 'HII_linelog.txt')
    # lime.save_line_log(sdss_spec.log, './sample_data/SDSS_flux_table.pdf')
    # lime.save_line_log(sdss_spec.log, './sample_data/SDSS_linelog.fits')
    # lime.save_line_log(sdss_spec.log, './sample_data/SDSS_linelog.xls')

    return ()

wave_spec, flux_spec = np.loadtxt('./spec_HI_Asier.dat', unpack=True, skiprows=1)

HII_flux(wave_spec, flux_spec)