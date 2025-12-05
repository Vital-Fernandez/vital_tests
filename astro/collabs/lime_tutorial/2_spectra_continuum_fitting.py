import lime

fname = '/home/vital/Astrodata/J0823_p2806.fits'
spatial_mask_file = './J0823_p2806_SN_mask.fits'
output_log = './J0823_p2806_lines_log.fits'

# Open the cube
cube = lime.Cube.from_file(fname, instrument='kcwi', redshift=0.04726, norm_flux=1e-6)
# cube.check.cube('H1_4861A', rest_frame=False)
#
# Central region spaxel
coords_HIIregion = (85, 85)
spec_HIIregion = cube.get_spectrum(coords_HIIregion[0], coords_HIIregion[1])
spec_HIIregion.fit.continuum(degree_list=[6,7,7,7], emis_threshold=[4,3,2, 1.5], smooth_scale=10, plot_steps=True)

# spec_HIIregion.plot.spectrum(show_err=True)
# spec_HIIregion.fit.bands('H1_4861A', err_from_bands=True, cont_source='adjacent')
# spec_HIIregion.plot.bands()

coords_Star = (27, 66)
spec_star = cube.get_spectrum(coords_Star[0], coords_Star[1])
spec_star.update_redshift(0)
spec_star.fit.continuum(degree_list=[6,7,7,7], emis_threshold=[4,3,2, 1.5], smooth_scale=30, plot_steps=True)

spec_star.fit.bands('H1_4861A', shape='abs', err_from_bands=True, cont_source='central')
spec_star.plot.bands(show_cont=True)
spec_star.fit.bands('H1_4861A', shape='abs', err_from_bands=True, cont_source='adjacent')
spec_star.plot.bands(show_cont=True)
spec_star.fit.bands('H1_4861A', shape='abs', err_from_bands=True, cont_source='fit')
spec_star.plot.bands(show_cont=True)


