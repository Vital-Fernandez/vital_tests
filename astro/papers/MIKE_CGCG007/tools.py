import numpy as np
from astropy.io import fits


def read_mike_fits(file_path):

    with fits.open(file_path) as hdul:
        wave_m, flux_m, err_m = hdul[1].data, hdul[2].data, hdul[3].data

    return wave_m, flux_m, err_m