from numpy import nan
import lime

fig_cfg = {'figure.dpi': 600,
           'figure.figsize': (16, 9),
           'axes.labelsize': 22,
           'legend.fontsize': 12,
           'legend.framealpha': 0.95,
           'xtick.labelsize': 20,
           'ytick.labelsize': 20}

lime.theme.set_style(fig_cfg=fig_cfg, colors_conf={'comp_width': 2})

# Declare the spectra and configuration files location:
fits_address = '/home/vital/PycharmProjects/ceers-data/data/spectra/CEERs_DR0.7/nirspec9/comb-mgrat/hlsp_ceers_jwst_nirspec_nirspec9-018294_comb-mgrat_v0.7_x1d-masked.fits'
conf_address = 'ceers_conf.toml'

# Fitting configuration
cfg = lime.load_cfg(conf_address)

# Create spectrum object (masking nan entries):
ceers18294 = lime.Spectrum.from_file(fits_address, 'nirspec', mask_flux_entries='nan', redshift=2.63645,
                                    units_flux='Jy')

# Plot the spectrum
# ceers18294.plot.spectrum()

# Change the units (and add a normalization)
ceers18294.unit_conversion(wave_units_out='AA', flux_units_out='FLAM', norm_flux=1e-22)

# Plot the
# ceers18294.plot.spectrum(rest_frame=False)

# # Getting the default line database in vacuum for the observed spectrum range
linesDF = lime.line_bands(wave_intvl=ceers18294, vacuum=True)
# print(linesDF)

# It seems there is a mask in the band of Hbeta and that was giving errors in the fitting... I need to check it out.
# But you can adjust the mask graphically with this command:
# Right click remove, middle click chane line notation, left click adjust the bands
# ceers18294.check.bands('ceers18294_bands_selection.txt', linesDF)

# Fit one line
ceers18294.fit.bands('H1_6565A_b', 'ceers18294_bands_selection.txt', cfg['default_line_fitting'])

# Plot the last fitting
ceers18294.plot.bands(fig_cfg=lime.theme.fig_defaults(), output_address='Halpha_ceers18294.png', y_scale='log')

# # Now we can measure the lines from our selection
# ceers6563.fit.frame('ceers6563_bands_selection.txt')
# ceers6563.plot.grid(output_address='ceers6563_line_grid.png')
# ceers6563.plot.spectrum(include_fits=True)
#
# # Save the measurements
# ceers6563.save_log('ceers1027_line_fluxes.txt')

