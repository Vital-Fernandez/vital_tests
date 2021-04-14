from astropy.io import fits

fits_addres = '/home/myfits.fits'

hdu_list = fits.open(fits_addres)
hdu_list[0] # primary hdu

hdu_list[0].header['DATE']
hdu_list[0].header[7]

hdr = hdu_list[0].header
                      # value           comment
hdr['targname'] = ('NGC121-a', 'the observation target')
hdr.comments['targname']

print(repr(hdr)) #use repr for a legible display

hdu_list.info() # Displays the information

hdu_list.close() # At the end

# Extension names "EXTNAME keywords" can be used to acces the fits data:
data = hdu_list['SCI'].data

# EXTVER can be used for extensions with the same name
data = hdu_list['SCI', 2].data

# Or use a context manager

with fits.open(fits_addres, lazy_load_hdus=False) as hdu_list:
    # lazy_load_hdus = False to make all the headers accesible after HDUList closes
    hdu_list.info()



