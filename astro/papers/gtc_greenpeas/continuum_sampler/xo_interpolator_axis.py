import numpy as np
import pyneb as pn
import scipy as spy
import exoplanet as xo

# Exoplanet 2D RegularGridInterpolator applied to PyNeb grids
print('\n -- Emissivity interpolation:')

# Declare ions
He1 = pn.RecAtom('He', 1)
H1 = pn.RecAtom('H', 1)

# Declare true values
Te_true, ne_true, wave_true = 14567.0, 275.0, 4026.0
emisTrue = He1.getEmissivity(Te_true, ne_true, wave=wave_true) / H1.getEmissivity(Te_true, ne_true, wave=4861)

# Define parameter grid
Te_range, ne_range = np.linspace(7000, 20000, 260), np.linspace(1, 1000, 100)
emisValues = He1.getEmissivity(Te_range, ne_range, wave=wave_true, product=True) / H1.getEmissivity(Te_range, ne_range, wave=4861, product=True)

# Define parameter 3D grid
linesNumber = 7
emisCube = np.ones((Te_range.size, ne_range.size, linesNumber))
for i in range(linesNumber):
    emisCube[:, :, i] = emisValues * i

# ---- Display results
print('1 True value', emisTrue)

# Scipy interpolation RegularGridInterpolator
spy_RGridInterp = spy.interpolate.RegularGridInterpolator((Te_range, ne_range), emisCube)
print('2 Scipy RegularGridInterpolator', spy_RGridInterp([[Te_true, ne_true]]))

# Exoplanet interpolation
coordB = np.stack(([Te_true], [ne_true]), axis=-1)
exop_interpAxis = xo.interp.RegularGridInterpolator([Te_range, ne_range], emisCube)
print('3 Exoplanet Axis Interpolation', exop_interpAxis.evaluate(coordB).eval())
