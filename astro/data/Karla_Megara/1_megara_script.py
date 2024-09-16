import urllib.request
from astropy.io import fits
from pathlib import Path
from astropy.wcs import WCS
from sys import stdout
import lime
import astropy


# Check version >= 0.9.96
print(lime.__version__)

# File location
path = Path('./')
cfgFile = path/'cfg 1.toml'
# obs_cfg = lime.load_cfg(cfgFile)
megara_cube_address = path/'NGC5471_datacube_LR-B_900_scale03_drp_nosky 1.fits'
bands_file_0 = path/'bands_a.txt'
output_lines_log_file = path/'NGC4571_lines.fits'
spatial_mask = path/'NGC5471_mask.fits'
output_lines_log_file = path/'NGC5471_lines.fits'

obs_cfg = lime.load_cfg(cfgFile)
redshift = obs_cfg['sample_data']['z']
norm_flux = obs_cfg['sample_data']['norm_flux']

ngc5471 = lime.Cube.from_file(megara_cube_address, instrument='megara', redshift=redshift)
ngc5471.unit_conversion(units_flux='FLAM', norm_flux=norm_flux)
# ngc5471.plot.cube(4341)

ngc5471.spatial_masking('H1_4861A', param='SN_line', contour_pctls=[93, 96, 99], output_address=spatial_mask)
# ngc5471.plot.cube('H1_4861A', masks_file=spatial_mask)
# ngc5471.check.cube('H1_4861A', masks_file=spatial_mask, rest_frame=True)

spaxel = ngc5471.get_spectrum(9, 9)
# spaxel.plot.spectrum(log_scale=False)
spaxel.fit.frame(bands=bands_file_0, fit_conf=cfgFile, line_detection=True, id_conf_prefix='MASK_0')
# spaxel.plot.spectrum(include_fits=True, rest_frame=True, log_scale=False)

