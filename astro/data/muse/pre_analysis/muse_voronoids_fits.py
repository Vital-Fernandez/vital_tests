import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, cm, colors
from astropy.wcs import WCS
from src.specsiser.print.plot import STANDARD_PLOT
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, image_array_binning
from astropy.io import fits


# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])
coordinates_keys_list = obsData['data_location']['wcs_key_list']

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']

# Plot set up
labelsDict = {'xlabel': r'RA',
              'ylabel': r'DEC'}
defaultConf = STANDARD_PLOT.copy()
defaultConf.update(labelsDict)
rcParams.update({})

bins_cords_fits_adress = fitsFolder/'VorbinCGCG'/'Bins_coordinates_map_CGCG007.fits'
bins_spec_fits_adress = fitsFolder/'VorbinCGCG'/'CGCG007_BinSpectra_linear.fits'

fits.info(bins_cords_fits_adress)
fits.info(bins_spec_fits_adress)

coord_hdr = fits.getheader(bins_cords_fits_adress, 'COORDINATES_MAP')
coord_data = fits.getdata(bins_cords_fits_adress, 'COORDINATES_MAP')

spec_hdr = fits.getheader(bins_spec_fits_adress, 'BIN_SPECTRA')
spec_data_flux = fits.getdata(bins_spec_fits_adress, 'BIN_SPECTRA')
spec_data_lamb = fits.getdata(bins_spec_fits_adress, 'LOGLAM')

flux_cube, err_cube = spec_data_flux['SPEC'], spec_data_flux['ESPEC']
wave = spec_data_lamb['LOGLAM']

binID = coord_data['BIN_ID'].astype(int)
binY, binX = coord_data['X_CUBE'].astype(int), coord_data['Y_CUBE'].astype(int)

print(binID)
