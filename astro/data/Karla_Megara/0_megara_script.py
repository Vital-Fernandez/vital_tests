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
path = Path('/home/vital/Astrodata/Karla_MEGARA/')
megara_cube_address = path/'NGC5471_datacube_LR-B_900_scale03_drp_nosky 1.fits'
bands_file = path/'bands_a.xlsx'
output_lines_log_file = path/'NGC4571_lines.fits'
spatial_mask = path/'NGC5471_mask.fits'

redshift = 0.00091
norm_flux = None

# Create LiMe cube
ngc5471 = lime.Cube.from_file(megara_cube_address, instrument='megara', redshift=redshift)
ngc5471.unit_conversion(flux_units_out='FLAM', norm_flux=norm_flux)
# ngc5471.check.cube('H1_4861A', rest_frame=True) #Hgamma
# ngc5471.plot.cube(4861) #Hbeta

# ngc5471.check.cube('H1_4861A', line_fg='H1_4861A', cont_pctls_fg=[93, 96, 99])
ngc5471.spatial_masking('H1_4861A', param='flux', contour_pctls=[93, 96, 99], output_address=spatial_mask)
ngc5471.check.cube('H1_4861A', masks_file=spatial_mask, rest_frame=True)

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