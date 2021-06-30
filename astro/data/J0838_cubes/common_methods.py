import numpy as np
import astropy.io.fits as astrofits

def import_fits_data(file_address, crval3, frame_idx=None):

    # Open fits file
    with astrofits.open(file_address) as hdul:
        data, hdr = hdul[frame_idx].data, hdul[frame_idx].header

    dw = hdr['CD3_3']
    w_min = crval3
    nPixels = hdr['NAXIS3']
    w_max = w_min + dw * nPixels
    wave = np.linspace(w_min, w_max, nPixels, endpoint=False)

    return wave, data, hdr