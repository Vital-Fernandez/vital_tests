import numpy as np
import pandas as pd
import src.specsyzer as sr
from pathlib import Path
from mpdaf.obj import Cube
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.stats import gaussian_fwhm_to_sigma
from mpdaf.obj import deg2sexa
from astropy.wcs import WCS
from photutils import detect_threshold, detect_sources
from astropy.convolution import Gaussian2DKernel
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.visualization import SqrtStretch

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

# Redshift correction
z = 0.0046
wave = wave / (1 + z)

# Get line flux region
lineLabel, lineWave = 'O3_5007A', 5007.0
lineRegions = np.array([lineWave-5, lineWave, lineWave+5])
lineIdcs = np.searchsorted(wave, lineRegions)
lineImage = cube[lineIdcs[0]:lineIdcs[2], :, :].sum(axis=0)
flux_image = lineImage.data.data

threshold = detect_threshold(flux_image, nsigma=5)

sigma = 3.0 * gaussian_fwhm_to_sigma  # FWHM = 3.
kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
kernel.normalize()
segm = detect_sources(flux_image, threshold, npixels=5, filter_kernel=kernel)
norm = ImageNormalize(stretch=SqrtStretch())

# Plot data
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12.5))
ax1.imshow(flux_image, origin='lower', cmap='Greys_r', norm=norm)
ax1.set_title('Data')
cmap = segm.make_cmap(random_state=12345)
ax2.imshow(segm, origin='lower', cmap=cmap)
ax2.set_title('Segmentation Image')
plt.show()

# # Plot line image map with coordinates
# labelsDict = {'xlabel': r'RA',
#               'ylabel': r'DEC',
#               'title': r'Galaxy CGCG007 $H\alpha$'}
# sr.plot.imageContour(flux_image, axes_conf=labelsDict)


