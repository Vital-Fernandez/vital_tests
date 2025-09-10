import numpy as np
import pyneb as pn
import scipy as spy
# import exoplanet as xo
import pandas as pd

# Declare ions
He1 = pn.RecAtom('He', 1)
H1 = pn.RecAtom('H', 1)

# Declare true values
Te_true, ne_true, wave_true, abund_true = 14567.0, 275.0, 4026.0, 7.45
emisTrue = abund_true * He1.getEmissivity(Te_true, ne_true, wave=wave_true) / H1.getEmissivity(Te_true, ne_true, wave=4861)
print('1 True value', emisTrue)

# Define parameter grid
Te_range, ne_range, abund_range = np.linspace(7000, 20000, 260), np.linspace(1, 1000, 100), np.arange(7.0, 8.0, 0.1)
emisValues = He1.getEmissivity(Te_range, ne_range, wave=wave_true, product=True) / H1.getEmissivity(Te_range, ne_range, wave=4861, product=True)

# Define parameter 3D grid
emisCube = np.ones((Te_range.size, ne_range.size, abund_range.size))
for i, abund_value in enumerate(abund_range):
    emisCube[:, :, i] = emisValues * abund_value

# Scipy interpolation interp2d
# print(emisValues.shape)
# spy_interp2d = spy.interpolate.interp2d(ne_range, Te_range, emisValues)
# print('2 Scipy interp2d', spy_interp2d(ne_range[:], Te_true).reshape(emisValues.shape))

# Scipy interpolation RegularGridInterpolator
spy_RGridInterp = spy.interpolate.RegularGridInterpolator((Te_range, ne_range, abund_range), emisCube)
print('3 Scipy RegularGridInterpolator', spy_RGridInterp([[Te_range, ne_range, abund_true]]))

# # Exoplanet interpolation
# exop_interp = xo.interp.RegularGridInterpolator([Te_range, ne_range, abund_range], emisCube, nout=1)
# coordB = np.stack(([Te_true], [ne_true], [abund_true]), axis=-1)
# print('4 Exoplanet interpolation', exop_interp.evaluate(coordB).eval())
#
# # --------------------------- Saving to a text file and loading grid ---------------------------------------------------
# output_DF = pd.Dataframe(columns=['Abund', 'Te', 'ne', 'O3_4363A'])
# output_file = 'D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/treatment/grid_file.txt'
#
# # output_DF['Abund'] =
# # output_DF['Te'] =
# output_DF['ne'] =
# output_DF['O3_4363A'] =
