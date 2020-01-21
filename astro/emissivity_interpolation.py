import numpy as np
import pyneb as pn
import scipy as spy
import exoplanet as xo

# Exoplanet 2D RegularGridInterpolator example
print('\n -- Exoplanet example interpolation')

x_range, y_range = np.linspace(-5, 5, 260), np.linspace(-1, 1, 100)
values = np.exp(-0.5 * (x_range ** 2)[:, None] - 0.5 * (y_range ** 2 / 0.5)[None, :] - x_range[:, None] * y_range[None, :])
interpolator = xo.interp.RegularGridInterpolator([x_range, y_range], values[:, :, None], nout=1)

i, j = 1, 0.5
print('True value', np.exp(-0.5 * (i ** 2) - 0.5 * (j ** 2 / 0.5) - i * j))

coordA = np.stack(([i], [j]), axis=-1)
print('Exoplanet interpolation', interpolator.evaluate(coordA).eval())

# Exoplanet 2D RegularGridInterpolator applied to PyNeb grids
print('\n -- Emissivity interpolation:')

# Declare ions
He1 = pn.RecAtom('He', 1)
H1 = pn.RecAtom('H', 1)

# Declare true values
Te_true, ne_true, wave_true = 14567.0, 275.0, 4026.0
emisTrue = He1.getEmissivity(Te_true, ne_true, wave=wave_true) / H1.getEmissivity(Te_true, ne_true, wave=4861)
print('True value', emisTrue)

# Define parameter grids
Te_range, ne_range = np.linspace(7000, 20000, 260), np.linspace(1, 1000, 100)
emisValues = He1.getEmissivity(Te_range, ne_range, wave=wave_true, product=True) / H1.getEmissivity(Te_range, ne_range, wave=4861, product=True)

# Scipy interpolation
spy_interp = spy.interpolate.interp2d(ne_range, Te_range, emisValues)
print('Scipy interpolation', spy_interp(ne_true, Te_true))

# Exoplanet interpolation
exop_interp = xo.interp.RegularGridInterpolator([Te_range, ne_range], emisValues[:, :, None], nout=1)
coordB = np.stack(([Te_true], [ne_true]), axis=-1)