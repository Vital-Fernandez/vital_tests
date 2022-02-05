import numpy as np
import pyneb as pn

# Compute the nSII density
T_low = 15000.0
ne_low = 250
S3 = pn.Atom('S', 3)
S2 = pn.Atom('S', 2)
O2 = pn.Atom('O', 2)
O3 = pn.Atom('O', 3)
mc_steps = 1000

# Declare observational data
objList = ['J1152', 'J0925'] # List of objects

# Fluxes of the emission lines normalized by Hbeta
data_dict = {'J1152': {'O2_3726_flux': 0.000,
                       'O2_3726_errF': 0.000,
                       'O2_3729_flux': 0.000,
                       'O2_3729_errF': 0.000},

             'J0925': {'O2_3726_flux': 0.000,
                       'O2_3726_errF': 0.000,
                       'O2_3729_flux': 0.000,
                       'O2_3729_errF': 0.000}} 

# Loop through the nights
for obj in objList:
    flux_dict = data_dict[obj]
    O2_3726A_dist = np.random.normal(flux_dict['O2_3726_flux'], flux_dict['O2_3726_errF'], mc_steps)
    O2_3729A_dist = np.random.normal(flux_dict['O2_3729_flux'], flux_dict['O2_3729_errF'], mc_steps)
    nOII = O2.getTemDen(O2_3726A_dist/O2_3729A_dist, tem=T_low, wave1=3726, wave2=3729)
    print(f'Galaxy {obj}: ne_OII = {nOII.mean():0.2f}+/-{nOII.std():0.2f} (cm^-3)')

