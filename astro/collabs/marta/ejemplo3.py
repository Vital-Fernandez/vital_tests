import numpy as np
from astropy.io import fits
import lime
import pandas as pd
from lime.io import _PARENT_BANDS

# Load spectrum
df = pd.read_csv('spec_9190-6102_1.txt')
wave = np.array(df.wave)
flux = np.array(df.flux)

# Galaxy redshift and the flux normalization
z_obj = 0
normFlux = 1

line = 'H1_6563A_b'
w_i = 6548
w_f = 6584
band_edges = np.array([w_i-60, w_i-25, w_i-10, w_f+10, w_f+25, w_f+60])

# Define a spectrum object
gp_spec = lime.Spectrum(wave, flux, redshift=z_obj, units_flux='1e20*FLAM')
# gp_spec.plot.spectrum(label='9190-6102')

# Run the fit
# gp_spec.fit.bands(line, band_edges)

# Plot the results from the last fitting
#gp_spec.plot.bands()

# Fit configuration           
# line = 'H1_6563A_b'
# fit_conf = {'H1_6563A_b': 'H1_6563A+N2_6584A+N2_6548A',
#             'N2_6548A_amp': {'expr': 'N2_6584A_amp/2.94'},
#             'N2_6548A_kinem': 'N2_6584A'}

# Second attempt including the fit configuration
gp_spec.fit.bands('H1_4861A_b', fit_conf='./fitting_conf.toml')
gp_spec.plot.bands()

gp_spec.fit.bands(line, band_edges, './fitting_conf.toml')
gp_spec.plot.bands()

# You can also save the fitting plot to a file
gp_spec.plot.spectrum()
