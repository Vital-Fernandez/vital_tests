import numpy as np
from astropy.io import fits
import lime

def wavelength_array_calculation(hdr):

    w_min = hdr['CRVAL1']
    dw = hdr['CDELT1']  # dw (Wavelength interval per pixel)
    pixels = hdr['NAXIS1']  # nw number of output pixels
    w_max = w_min + dw * pixels
    wave = np.linspace(w_min, w_max, pixels, endpoint=False)

    return wave

file1 = '/home/vital/Downloads/ADP.2014-05-15T19_19_57.533.fits'
file2 = '/home/vital/Downloads/ADP.2014-05-15T19_19_57.507.fits'

# Checking the format of the file 1
print(fits.info(file1))

# Checking the format of the file 2
print(fits.info(file2))

# with fits.open(file1) as hdul:
#     data1, header1 = hdul[0].data, hdul[0].header
#     err1 = hdul[1].data
# #
with fits.open(file1) as hdul:
    data2, header2 = hdul[1].data, hdul[1].header

# wave1 = wavelength_array_calculation(header1)
# spec1 = lime.Spectrum(wave1, data1[0], redshift=0.283241)
# spec1.plot_spectrum(frame='rest')
#
# wave2 = wavelength_array_calculation(header2)
# spec2 = lime.Spectrum(wave2, data2[1], redshift=0.283241)
# spec2.plot_spectrum(frame='rest')

with fits.open(file2) as hdul:
    data2, header2 = hdul[1].data, hdul[1].header

wave = data2['WAVE']*10
flux = data2['FLUX'][0]
3