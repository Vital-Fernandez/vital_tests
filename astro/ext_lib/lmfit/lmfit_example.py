import matplotlib.pyplot as plt
import numpy as np
from lmfit.models import ExponentialModel, GaussianModel

dat = np.loadtxt('NIST_Gauss2.dat')
x = dat[:, 1]
y = dat[:, 0]

model = ExponentialModel(prefix='exp_')
model.set_param_hint('amplitude', value=10)
model.set_param_hint('decay', value=10)

model += GaussianModel(prefix='g1_')
model.set_param_hint('g1_center', value=105, min=75, max=125)
model.set_param_hint('g1_sigma', value=15, min=3)
model.set_param_hint('g1_amplitude', value=2000, min=10)

model += GaussianModel(prefix='g2_')
model.set_param_hint('g2_center', value=155, min=125, max=175)
model.set_param_hint('g2_delta_sigma', value=1.5, min=0.8)
model.set_param_hint('g2_sigma', expr='g2_delta_sigma*g1_sigma')
model.set_param_hint('g2_amplitude', value=2000, min=10)

pars = model.make_params()

init = model.eval(pars, x=x)
out = model.fit(y, pars, x=x)

print(out.fit_report(min_correl=0.5))

fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.8))
axes[0].plot(x, y, 'b')
axes[0].plot(x, init, 'k--', label='initial fit')
axes[0].plot(x, out.best_fit, 'r-', label='best fit')
axes[0].legend(loc='best')

comps = out.eval_components(x=x)
axes[1].plot(x, y, 'b')
axes[1].plot(x, comps['g1_'], 'g--', label='Gaussian component 1')
axes[1].plot(x, comps['g2_'], 'm--', label='Gaussian component 2')
axes[1].plot(x, comps['exp_'], 'k--', label='Exponential component')
axes[1].legend(loc='best')

plt.show()