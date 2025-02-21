from astroquery.sdss import SDSS
from astropy import coordinates as coords
pos = coords.SkyCoord('09h58m38.4196s +20d25m7.63s', frame='icrs')
xid = SDSS.query_region(pos, radius='5 arcsec', spectro=True)
sp = SDSS.get_spectra(matches=xid)
im = SDSS.get_images(matches=xid, band='g')


from pathlib import Path

import numpy as np
import pandas as pd
import lime

# from astropy.io import fits
#
# fits_name = '/home/vital/Astrodata/JWST-GO-05554/MAST_2025-02-19T22_09_23.544Z/MAST_2025-02-19T22_09_23.544Z/JWST/jw02424-o902_t002_miri/jw02424-o902_t002_miri_ch4-shortmediumlong_s3d.fits'
# hdr = fits.getheader(fits_name, ext=0)
# hdr['TARGNAME']
# for key, value in hdr.items():
#     print(key, value)

# # Load the data configuration
# cfg = lime.load_cfg('cfg.toml')
#
# # Declare the data
# data_folder = Path(cfg['meta']['data_folder'])
#
#
# cube = lime.Cube.from_file()
# print(data_folder)