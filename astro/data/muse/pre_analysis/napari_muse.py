import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from mpdaf.obj import Cube
# import napari


# Load data
files_list = ['CGCG007.fits', 'UGC5205.fits']
files_path = Path('D:/Google drive/Astrophysics/Datos/MUSE - Amorin')
files_address = list(files_path/file for file in files_list)
cube = Cube(filename=str(files_address[0]))

# Cube shape (lambda, Y, X) // Header shape (1, 2, 3) = (RA---TAN, DEC--TAN, AWAV) = (X, Y, lambda) // QfitsView (X, Y)
cube_size = cube.shape

# Fits properties
cube.info()

# Fits header
hdr = cube.data_header

# Reconstruct wave:
cube.wave.info()
dw = hdr['CD3_3']
w_min = hdr['CRVAL3']
nPixels = hdr['NAXIS3']
w_max = w_min + dw * nPixels
wave = np.linspace(w_min, w_max, nPixels, endpoint=False)

# with napari.gui_qt():
#     viewer = napari.view_image(cube.data)

