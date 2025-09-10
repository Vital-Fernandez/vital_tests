from numpy import nan
import lime

# Declare the spectra and configuration files location:
fits_address = 'hlsp_ceers_jwst_nirspec_nirspec10-001027_comb-mgrat_v0.7_x1d-masked.fits'
conf_address = 'ceers_conf.toml'
cfg = lime.load_cfg(conf_address)

# Create spectrum object (masking nan entries):
ceers1027 = lime.Spectrum.from_file(fits_address, instrument='nirspec', mask_flux_entries=[nan], redshift=7.8189)

# Plot the spectrum
ceers1027.plot.spectrum()

# Change the units (and add a normalization)
ceers1027.unit_conversion(units_wave='A', units_flux='Flam', norm_flux=1e-22)

# Plot the
ceers1027.plot.spectrum(rest_frame=True)

# Measuring a line from the database and show fit (last one made)
ceers1027.fit.bands(4861)
ceers1027.plot.bands()
ceers1027.plot.bands(rest_frame=True)

# # Getting the default line database in vacuum for the observed spectrum range
linesDF = lime.line_bands(wave_intvl=ceers1027, vacuum=True)
print(linesDF)
#
# # # Fitting several lines from this database using the LiMe notation
# ceers1027.fit.frame(linesDF, line_list=['H1_4863A', 'O3_4960A', 'O3_5008A'])
#
# # Plot the fittings as a grid and in the spectrum (we have two O3_5007A from before)
# ceers1027.plot.grid()
# ceers1027.plot.spectrum(include_fits=True)

# It seems there is a mask in the band of Hbeta and that was giving errors in the fitting... I need to check it out.
# But you can adjust the mask graphically with this command:
# Right click remove, middle click chane line notation, left click adjust the bands
ceers1027.check.bands('ceers1027_bands_selection.txt', linesDF)

# Clearing the log just in case
ceers1027.log = ceers1027.log.iloc[0:0]

# Now we can measure the lines from our selection
ceers1027.fit.frame('ceers1027_bands_selection.txt', fit_conf=cfg)
ceers1027.plot.grid(output_address='ceers1027_line_grid.png')
ceers1027.plot.spectrum(include_fits=True)

# Save the measurements
ceers1027.save_log('ceers1027_line_fluxes.txt')
