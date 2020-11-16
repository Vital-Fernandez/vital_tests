import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp, loadtxt, pi, sqrt
from lmfit import Model
from lmfit.models import GaussianModel, LinearModel, PolynomialModel, height_expr

data = loadtxt('C:/Users/Vital/OneDrive/Desktop/model1d_gauss.dat')
x = data[:, 0]
y = data[:, 1] + 0.25 * x - 1.0
SQRT2PI = np.sqrt(2*pi)

def gauss_area(sigma_true, amp_true):
    # return np.sqrt(2 * np.pi * sigma_true ** 2) * amp_true
    return amp_true * SQRT2PI * sigma_true

def line(x, slope, intercept):
    """a line"""
    return slope * x + intercept

def gaussian_height(x, amp, cen, wid):
    """1-d gaussian curve : gaussian(x, amp, cen, wid)"""
    return amp * np.exp(-0.5 * (((x-cen)/wid) * ((x-cen)/wid)))

def gaussian_a(x, amp, cen, wid):
    """1-d gaussian: gaussian(x, amp, cen, wid)"""
    return (amp / (SQRT2PI * wid)) * exp(-(x-cen)**2 / (2*wid**2))


# def gaussian(x, amp, cen, wid):
#     """1-d gaussian: gaussian(x, amp, cen, wid)"""
#     return (amp / (sqrt(2*pi) * wid)) * exp(-(x-cen)**2 / (2*wid**2))

# # ---------------- USING YOUR OWN MODELS ORIGINAL
# mod = Model(gaussian_a) + Model(line)
# pars = mod.make_params(amp=5, cen=5, wid=1, slope=0, intercept=1)
# pars.pretty_print()
# result = mod.fit(y, pars, x=x)

# ---------------- USING YOUR OWN MODELS PLAYING
idx_max = np.argmax(y)
mod_h = Model(line)
mod_h.prefix = 'myCont_'
mod_h += Model(gaussian_height, prefix='gauss_')
mod_h.set_param_hint('gauss_amp', value=y[idx_max], min=0)
mod_h.set_param_hint('gauss_cen', value=x[idx_max], putoNewwile=True)
mod_h.set_param_hint('gauss_wid', value=1, min=0)
mod_h.set_param_hint('myCont_intercept', value=0)
mod_h.set_param_hint('myCont_slope', value=1)
pars = mod_h.make_params()
pars.pretty_print()
result = mod_h.fit(y, pars, x=x)
print(gauss_area(result.params['gauss_wid'].value, result.params['gauss_amp'].value))



# mod_a = Model(gaussian_a, prefix='gauss_') + Model(line, prefix='cont_')
# mod_a.set_param_hint('gauss_amp', value=5, min=0)
# mod_a.set_param_hint('gauss_cen', value=5, putoNewwile=True)
# mod_a.set_param_hint('gauss_wid', value=1, min=0)
# mod_a.set_param_hint('cont_intercept', value=0)
# mod_a.set_param_hint('cont_slope', value=1)
# pars = mod_a.make_params()
# pars.pretty_print()
# result = mod_a.fit(y, pars, x=x)


# ---------------- USING DEFAULT MODELS
# continuumModel = LinearModel(prefix='cont_')
# narrowModel = GaussianModel(prefix=f'gauss_')
# mod = continuumModel + narrowModel
# mod.set_param_hint('gauss_amplitude', value=5, min=0, coso=3)
# mod.set_param_hint('gauss_center', value=5, min=0, coso=3)
# mod.set_param_hint('gauss_sigma', value=1, min=0, coso=3)
# mod.set_param_hint('cont_intercept', value=0)
# mod.set_param_hint('cont_slope', value=1)
# pars = mod.make_params()
# result = mod.fit(y, pars, x=x)
# print(1/SQRT2PI)
# '0.3989423*gauss_amplitude/max(2.220446049250313e-16, gauss_sigma)'
# np.sqrt(2 * np.pi * sigmaTrue ** 2) * ampTrue
# ---------------- RUN THE MODEL

# ---------------- DISPLAY THE RESULTS
print(result.fit_report())

plt.plot(x, y, 'bo')
plt.plot(x, result.init_fit, 'k--', label='initial fit')
plt.plot(x, result.best_fit, 'r-', label='best fit')
plt.legend(loc='best')
plt.show()



# """
# Fit Using Inequality Constraint
# ===============================
#
# Sometimes specifying boundaries using ``min`` and ``max`` are not sufficient,
# and more complicated (inequality) constraints are needed. In the example below
# the center of the Lorentzian peak is constrained to be between 0-5 away from
# the center of the Gaussian peak.
#
# See also: https://lmfit.github.io/lmfit-py/constraints.html#using-inequality-constraints
# """
#
# import matplotlib.pyplot as plt
# import numpy as np
#
# from lmfit import Minimizer, Parameters, report_fit
# from lmfit.lineshapes import gaussian, lorentzian
#
#
# def residual(pars, x, data):
#     model = (gaussian(x, pars['amp_g'], pars['cen_g'], pars['wid_g']) +
#              lorentzian(x, pars['amp_l'], pars['cen_l'], pars['wid_l']))
#     return model - data
#
#
# ###############################################################################
# # Generate the simulated data using a Gaussian and Lorentzian line shape:
# np.random.seed(0)
# x = np.linspace(0, 20.0, 601)
#
# data = (gaussian(x, 21, 6.1, 1.2) + lorentzian(x, 10, 9.6, 1.3) +
#         np.random.normal(scale=0.1, size=x.size))
#
# ###############################################################################
# # Create the fitting parameters and set an inequality constraint for ``cen_l``.
# # First, we add a new fitting  parameter ``peak_split``, which can take values
# # between 0 and 5. Afterwards, we constrain the value for ``cen_l`` using the
# # expression to be ``'peak_split+cen_g'``:
#
# pfit = Parameters()
# pfit.add(name='amp_g', value=10)
# pfit.add(name='amp_l', value=10)
# pfit.add(name='cen_g', value=5)
# pfit.add(name='peak_split', value=2.5, min=0, max=5, vary=True)
# pfit.add(name='cen_l', expr='peak_split+cen_g')
# pfit.add(name='wid_g', value=1)
# pfit.add(name='wid_l', expr='wid_g')
#
# mini = Minimizer(residual, pfit, fcn_args=(x, data))
# out = mini.leastsq()
# best_fit = data + out.residual
#
# ###############################################################################
# # Performing a fit, here using the ``leastsq`` algorithm, gives the following
# # fitting results:
# report_fit(out.params)
#
# ###############################################################################
# # and figure:
# plt.plot(x, data, 'bo')
# plt.plot(x, best_fit, 'r--', label='best fit')
# plt.legend(loc='best')
# plt.show()
