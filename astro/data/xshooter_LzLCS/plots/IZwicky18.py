from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import src.specsiser as sr
import pyneb as pn


# spec_address = '/home/vital/Astro-data/Observations/IZW18_Blue_cr_f_t_w_e__test1_fglobal.fits'
#
# wave, data, header = sr.import_fits_data(spec_address, instrument='ISIS', frame_idx=0)
# norm_flux = 1e-17
# z_obj = 0.0
#
# lm = sr.LineMesurer(wave, data[0], normFlux=norm_flux, redshift=z_obj)
# lm.plot_spectrum(specLabel='IZwicky18 Blue arm ISIS spectrum')

spec_address = '/home/vital/Dropbox/Astrophysics/Seminars/LzLCS/spec-0266-51630-0100.fits'

wave, data, header = sr.import_fits_data(spec_address, instrument='SDSS', frame_idx=0)
norm_flux = 1e-17
z_obj = 0.0

lm = sr.LineMesurer(wave, data['flux']*norm_flux, normFlux=norm_flux, redshift=z_obj)
lm.plot_spectrum(specLabel='CGCG007-025 SLOAN spectrum')