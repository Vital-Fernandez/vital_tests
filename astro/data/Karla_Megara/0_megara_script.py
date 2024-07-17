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
# cfgFile = path+'cfg 1.toml'
# obs_cfg = lime.load_cfg(cfgFile)
megara_cube_address = path/'NGC5471_datacube_LR-B_900_scale03_drp_nosky.fits'
bands_file = path/'bands_a.xlsx'
output_lines_log_file = path/'NGC4571_lines.fits'
spatial_mask = path/'NGC5471_mask.fits'

redshift = 0.00091
norm_flux = None

# Create LiMe cube
ngc5471 = lime.Cube.from_file(megara_cube_address, instrument='megara', redshift=redshift)
ngc5471.unit_conversion(units_flux='FLAM', norm_flux=norm_flux)
#ngc5471.plot.cube(4341) #Hgamma
# ngc5471.plot.cube(4861) #Hbeta

#Plot the emission of Hbeta versus other line ratios in the same CUBE
# ngc5471.check.cube(4862, line_fg=5007, cont_pctls_fg=[85, 90, 95, 99])
# ngc5471.plot.cube(5007, line_fg=4959)
# ngc5471.spatial_masking('O3_4363A', param='SN_line', contour_pctls=[93, 96, 99], output_address=spatial_mask)
ngc5471.check.cube('H1_4861A', masks_file=spatial_mask, rest_frame = True)
ngc5471.check.cube('H1_4861A', lines_file=output_lines_log_file, rest_frame=True)
#Please continue with the Tutorial 4 and 5 of LIME.
#The goal is to make some emission lines: Hbeta and [O III] 5007. It possible to derive make a map
#with the extiction and electron temperature, E(b-V) and Te[O III] for this CUBE.
#Please be sure to use the most updated version of LIME.
# spaxel = ngc5471.get_spectrum(9, 9)
# spaxel.plot.spectrum(log_scale=False)
# spaxel.fit.frame(bands_file, cfgFile, line_detection=True, id_conf_prefix = 'MASK_0')
# spaxel.plot.spectrum(include_fits=True, rest_frame=True, log_scale=True)
# spaxel.load_frame(output_lines_log_file, page='9-9_LINELOG')
# spaxel.plot.grid()