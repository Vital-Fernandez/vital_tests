from astropy.io import fits
import lime

file1 = './SDSS020356_VIS.fits'
mask_file = './xshooter_mask_0356.txt'
cfg_file = './config.txt'

data = fits.getdata(file1, ext=1)
hdr = fits.getheader(file1, ext=1)
obsData = lime.load_cfg(cfg_file)
mask = lime.load_lines_log(mask_file)

wave = data['WAVE'][0]*10
flux = data['FLUX'][0]

norm_obj = 1e-20
z_obj = 0.188

print(mask.loc['O3_4959A'].values)

gp_spec = lime.Spectrum(wave, flux, redshift=z_obj, norm_flux=norm_obj, crop_waves=(5440, 9900))
gp_spec.plot_spectrum(frame='rest')

gp_spec.fit_from_wavelengths('O3_4959A', mask.loc['O3_4959A'].values)
gp_spec.display_results(fit_report=True)

gp_spec.fit_from_wavelengths('O3_4959A_b', mask.loc['O3_4959A'].values, obsData['two_component_line_fitting'])
gp_spec.display_results(fit_report=True)

