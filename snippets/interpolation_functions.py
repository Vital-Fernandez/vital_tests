import numpy as np
import pyneb_diagnostic_diagram as pn
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
print('1 True value', emisTrue)

# Define parameter grid
Te_range, ne_range = np.linspace(7000, 20000, 260), np.linspace(1, 1000, 100)
emisValues = He1.getEmissivity(Te_range, ne_range, wave=wave_true, product=True) / H1.getEmissivity(Te_range, ne_range, wave=4861, product=True)

# Define parameter 3D grid
linesNumber = 7
emisCube = np.ones((Te_range.size, ne_range.size, linesNumber))
for i in range(linesNumber):
    emisCube[:, :, i] = emisValues * i

# Scipy interpolation interp2d
spy_interp2d = spy.interpolate.interp2d(ne_range, Te_range, emisValues)
print('2 Scipy interp2d', spy_interp2d(ne_true, Te_true))

# Scipy interpolation RegularGridInterpolator
spy_RGridInterp = spy.interpolate.RegularGridInterpolator((Te_range, ne_range), emisCube)
print('3 Scipy RegularGridInterpolator', spy_RGridInterp([[Te_true, ne_true]]))

# Exoplanet interpolation
exop_interp = xo.interp.RegularGridInterpolator([Te_range, ne_range], emisValues[:, :, None], nout=1)
coordB = np.stack(([Te_true], [ne_true]), axis=-1)
print('4 Exoplanet interpolation', exop_interp.evaluate(coordB).eval())

# Exoplanet interpolation
exop_interpAxis = xo.interp.RegularGridInterpolator([Te_range, ne_range], emisCube.reshape((Te_range.size, ne_range.size, -1)))
print('5 Exoplanet Axis Interpolation', exop_interpAxis.evaluate(coordB).eval())


# import numpy as np
# import pyneb as pn
# import exoplanet as xo
# import theano.tensor as T
#
# # Exoplanet 2D RegularGridInterpolator applied to PyNeb grids
# print('\n -- Emissivity interpolation:')
#
# # Declare ions
# He1 = pn.RecAtom('He', 1)
# H1 = pn.RecAtom('H', 1)
#
# # Declare true values
# Te_true, ne_true, wave_true = 14567.0, 275.0, 4026.0
# emisTrue = He1.getEmissivity(Te_true, ne_true, wave=wave_true) / H1.getEmissivity(Te_true, ne_true, wave=4861)
# print('1 True value', emisTrue)
#
# # Define parameter grid
# Te_range, ne_range = np.linspace(7000, 20000, 260), np.linspace(1, 1000, 100)
# emisValues = He1.getEmissivity(Te_range, ne_range, wave=wave_true, product=True) / H1.getEmissivity(Te_range, ne_range, wave=4861, product=True)
#
# # Exoplanet interpolation
# exop_interp = xo.interp.RegularGridInterpolator([Te_range, ne_range], emisValues[:, :, None], nout=1)
# coordB = np.stack(([Te_true], [ne_true]), axis=-1)
# coordB = np.stack(([Te_true], [ne_true]), axis=-1)
# print('4 Exoplanet interpolation', exop_interp.evaluate(coordB).eval())
#
# Te, ne = T.dscalar('Te'), T.dscalar('ne')
# coord = T.stack([Te, ne], axis=-1)
# func = exop_interp.evaluate(coord)
# func.eval()






