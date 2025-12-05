import lime

# Change the database transtitions to vacuum
lime.lineDB.set_database(vacuum_waves=True)

# State the data location
fname = '/home/vital/Downloads/bluejay-south-v3_g235m-f170lp_1810_10314.spec.fits'
cfg_file = 'bluejay_cfg.toml'
obj_bands_path = 'bluejay-south_bands.txt'

# Open the file
spec = lime.Spectrum.from_file(fname, instrument='nirspec_grizli', redshift=2.09987)
spec.unit_conversion('AA', 'FLAM')
# spec.plot.spectrum(rest_frame=True)

# Generate the line bands (both obj_cfg_prefix give same results)
obj_bands = spec.retrieve.lines_frame(band_vsigma=140, fit_cfg=cfg_file, obj_cfg_prefix='Halpha_plus_nitrogen')
lime.save_frame(obj_bands_path, obj_bands)
# spec.plot.spectrum(bands=obj_bands, rest_frame=True)

# Fit the line
spec.fit.bands('H1_6563A_b', obj_bands_path, cfg_file, obj_cfg_prefix='Halpha_plus_nitrogen')
spec.plot.bands()

