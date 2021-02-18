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

# Comparison
print('Exoplanet interpolation', exop_interp.evaluate(coordB).eval())

# import numpy as np
# import matplotlib.pyplot as plt
# import pymc3 as pm
# import exoplanet as xo
# import theano.tensor as tt
#
# # np.random.seed(42)
# # x = np.sort(np.random.uniform(-5, 5, 25))
# # points = [x]
# # values = x ** 3 - x ** 2
# #
# # interpolator = xo.interp.RegularGridInterpolator(points, values[:, None])
# #
# # t = np.linspace(-6, 6, 5000)
# # plt.plot(t, interpolator.evaluate(t[:, None]).eval(), label="interpolation")
# # plt.plot(x, values, "o", label="control points")
# # plt.xlabel("x")
# # plt.ylabel("y")
# # plt.legend(fontsize=12)
# # plt.show()
#
# # truth = 45.0
# # data_sd = 8.0
# # data_mu = truth + data_sd * np.random.randn()
# #
# # with pm.Model() as model:
# #
# #     # The value passed into the interpolator must have the shape
# #     # (ntest, ndim), but in our case that is (1, 1)
# #     xval = pm.Uniform("x", lower=-8, upper=8, shape=(1, 1))
# #
# #     # Evaluate the interpolated model and extract the scalar value
# #     # we want
# #     mod = pm.Deterministic("y", interpolator.evaluate(xval)[0, 0])
# #
# #     # The usual likelihood
# #     pm.Normal("obs", mu=mod, sd=data_sd, observed=data_mu)
# #
# #     # Sampling!
# #     trace = pm.sample(draws=1000, tune=2000, step_kwargs=dict(target_accept=0.9), cores=1)
# #
# # t = np.linspace(-6, 6, 5000)
# # plt.plot(trace["x"][:, 0, 0], trace["y"], ".k", alpha=0.1, label="posterior samples")
# # plt.axhline(truth, color="k", lw=3, alpha=0.5, label="truth")
# # plt.plot(t, interpolator.evaluate(t[:, None]).eval(), label="interpolation")
# # plt.xlabel("x")
# # plt.ylabel("y")
# # plt.legend(fontsize=12)
# # plt.show()
# # print('Hi')
#
# points = np.array([np.linspace(-5, 5, 50), np.linspace(-1, 1, 25)])
# values = np.exp(-0.5 * (points[0] ** 2)[:, None] - 0.5 * (points[1] ** 2 / 0.5)[None, :] - points[0][:, None] * points[1][None, :])
#
# interpolator = xo.interp.RegularGridInterpolator(points, values[:, :, None], nout=1)
#
# i, j = 1, 0.5
# k = np.exp(-0.5 * (i ** 2) - 0.5 * (j ** 2 / 0.5) - i * j)
# cord = np.array([i, j])
#
# print('This value', interpolator.evaluate(cord[:, :, None]).eval())
#
# data_mu = 0.6
# data_sd = 0.1
#
# cord = [i, j]
#
# with pm.Model() as model:
#
#     xval = pm.Uniform("x", lower=-5, upper=5, shape=(1,))
#     yval = pm.Uniform("y", lower=-1, upper=1, shape=(1,))
#     xtest = tt.stack([xval, yval], axis=-1)
#
#     mod = interpolator.evaluate(xtest)
#
#     # The usual likelihood
#     pm.Normal("obs", mu=mod, sd=data_sd, observed=data_mu)
#
#     # Sampling!
#     trace = pm.sample(draws=4000, tune=4000, step_kwargs=dict(target_accept=0.9), chains=2, cores=1)