# H1_6563A_w1_sigma = '>1.5*H1_6563A_sigma'
#
# param_ref = '>1.5*H1_6563A_sigma'
#
# # Create additional parameter
# ineq_name = f'{param_ref}_ineq'
# ineq_conf = dict(value=param_ref)
#
# # Prepare definition of param:
# new_expresion = param_ref.replace(param_ref, ineq_name)
# new_expresion = new_expresion.replace('<', '').replace('>', '')
# print(new_expresion)
#
# param_ref = '>1.5*H1_6563A_sigma'

import matplotlib.pyplot as plt
import numpy as np

from lmfit import Minimizer, Parameters, report_fit
from lmfit.lineshapes import gaussian, lorentzian


def residual(pars, x, data):
    model = (gaussian(x, pars['amp_g'], pars['cen_g'], pars['wid_g']) +
             lorentzian(x, pars['amp_l'], pars['cen_l'], pars['wid_l']))
    return model - data

np.random.seed(0)
x = np.linspace(0, 20.0, 601)

data = (gaussian(x, 21, 6.1, 1.0) + lorentzian(x, 10, 9.6, 1.3) +
        np.random.normal(scale=0.1, size=x.size))

pfit = Parameters()
# pfit.add(name='amp_g', value=10)
# pfit.add(name='amp_l', value=10)
# pfit.add(name='cen_g', value=5)
# pfit.add(name='peak_split', value=2.5, min=0, max=5, vary=True)
# pfit.add(name='cen_l', expr='peak_split+cen_g')
# pfit.add(name='wid_g', value=1)
# pfit.add(name='wid_l', expr='wid_g')

pfit.add(name='amp_g', value=10)
pfit.add(name='amp_l', value=10)
pfit.add(name='cen_g', value=5)
pfit.add(name='cen_l', value=10)
pfit.add(name='wid_g', value=1)
pfit.add(name='wid_l_ineq', value=0.5, min=1)
pfit.add(name='wid_l', expr='wid_l_ineq*wid_g')

mini = Minimizer(residual, pfit, fcn_args=(x, data))
out = mini.leastsq()
best_fit = data + out.residual

report_fit(out.params)

plt.plot(x, data, 'bo')
plt.plot(x, best_fit, 'r--', label='best fit')
plt.legend(loc='best')
plt.show()