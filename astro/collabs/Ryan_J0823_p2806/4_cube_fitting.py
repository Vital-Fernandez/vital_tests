import lime

# Input files
fname = '/home/vital/Astrodata/J0823_p2806_v2.fits'
cfgname = './example_cfg.toml'
combine_mask_file = './J0823_p2806_SN_mask.fits'

# Load the cube
cube = lime.Cube.from_file(fname, instrument='kcwi', redshift=0.04726)

# # Fit the lines on the spatial masks
# output_log = './J0823_p2806_FG_O2_60p_single_comp_log.fits'
# cube.fit.spatial_mask(combine_mask_file, output_log, fit_cfg=cfgname, mask_list=['FG_O2_60p'], line_detection=True,
#                       err_from_bands=True, cont_source='adjacent')

output_log = './J0823_p2806_Temp_O3_99p_log.fits'
cube.fit.spatial_mask(combine_mask_file, output_log, fit_cfg=cfgname, mask_list=['Temp_O3_99p'], line_detection=True,
                      err_from_bands=True, cont_source='adjacent')

# Check the results
cube.check.cube('O2_3726A', masks_file=combine_mask_file, min_pctl_bg=55, rest_frame=True, masks_cmap='viridis',
                fname=output_log)
