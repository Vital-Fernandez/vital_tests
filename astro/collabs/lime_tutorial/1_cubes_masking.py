import lime

fname = '/home/vital/Astrodata/J0823_p2806.fits'
spatial_O2nebular_mask_file = './J0823_p2806_SN_O2_3726A_mask.fits'
spatial_O3nebular_mask_file = './J0823_p2806_SN_O3_5007A_mask.fits'
spatial_O3auroral_mask_file = './J0823_p2806_SN_O3_4363A_mask.fits'
spatial_He2_mask_file = './J0823_p2806_SN_He2_4686A_mask.fits'

# Open the cube
cube = lime.Cube.from_file(fname, instrument='kcwi', redshift=0.04726)

cube.spatial_masking('O2_3726A', param='SN_cont', contour_pctls=[55, 60, 75, 90], fname=spatial_O2nebular_mask_file)
cube.check.cube('O2_3726A', masks_file=spatial_O2nebular_mask_file, min_pctl_bg=55, rest_frame=True, masks_cmap='viridis')

cube.spatial_masking('O2_3726A', param='SN_line', contour_pctls=[55, 60], fname=spatial_O2nebular_mask_file)
cube.check.cube('O2_3726A', masks_file=spatial_O2nebular_mask_file, min_pctl_bg=55, rest_frame=True, masks_cmap='viridis')

cube.spatial_masking('O3_4363A', param='SN_line', contour_pctls=[98, 98.5, 99.5], fname=spatial_O3auroral_mask_file)
cube.check.cube('O3_5007A', masks_file=spatial_O3auroral_mask_file, min_pctl_bg=60, rest_frame=True, masks_cmap='viridis')
