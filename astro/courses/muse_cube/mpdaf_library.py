import numpy as np
from mpdaf.obj import Cube

def import_fits_data(file_address):

    cube = Cube(filename=str(file_address))
    header = cube.data_header

    cube.wave.info()
    dw = header['CD3_3']
    w_min = header['CRVAL3']
    nPixels = header['NAXIS3']
    w_max = w_min + dw * nPixels
    wave = np.linspace(w_min, w_max, nPixels, endpoint=False)

    return wave, cube, header