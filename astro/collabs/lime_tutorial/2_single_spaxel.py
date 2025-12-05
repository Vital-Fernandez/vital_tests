import lime

# Input files
fname = '/home/vital/J0823_p2806.fits'
cfgname = './example_cfg.toml'
kcwi_fname = './kcwi_bands.txt'
spatial_O3auroral_mask_file = './J0823_p2806_SN_O3_4363A_mask.fits'

# Output files
sp_mask0_bands = './sp_mask0_bands.txt'
sp_mask1_bands = './sp_mask1_bands.txt'
output_log = './J0823_p2806_lines_log.fits'

# Load the cube
cube = lime.Cube.from_file(fname, instrument='kcwi', redshift=0.04726)

# Extract the spectrum from a high ionization region
spec = cube.get_spectrum(85, 85)

# Check the bands
spec.check.bands('./kcwi_fname.txt', show_continua=True)

# Generate the bands for each mask
spec_bands0 = spec.retrieve.lines_frame(band_vsigma=70, fit_cfg=cfgname, obj_cfg_prefix='MASK_0', ref_bands=kcwi_fname)
spec_bands1 = spec.retrieve.lines_frame(band_vsigma=70, fit_cfg=cfgname, obj_cfg_prefix='MASK_1', ref_bands=kcwi_fname)
lime.save_frame(sp_mask0_bands, spec_bands0)
lime.save_frame(sp_mask1_bands, spec_bands1)

# Fit the object continuum
spec.fit.continuum(degree_list=[5,6,7], emis_threshold=[4, 3, 2], smooth_scale=5, plot_steps=False)

# Matched the peaks observed with the bands
matched_bands = spec.infer.peaks_troughs(spec_bands0, sigma_threshold=3, plot_steps=False)
spec.plot.spectrum(bands=matched_bands, show_cont=True)

# Fit the lines
spec.fit.frame(matched_bands, fit_cfg=cfgname, obj_cfg_prefix=None, err_from_bands=True, cont_source='adjacent')
spec.plot.grid()


cube.fit.spatial_mask(spatial_O3auroral_mask_file, output_log, fit_cfg=cfgname, mask_list=['MASK_0', 'MASK_1'],
                      line_detection=True, err_from_bands=True, cont_source='adjacent')

cube.check.cube('O3_5007A', masks_file=spatial_O3auroral_mask_file, min_pctl_bg=55, rest_frame=True, masks_cmap='viridis',
                fname=output_log)
