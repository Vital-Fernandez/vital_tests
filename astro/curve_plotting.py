import numpy as np
import pymc3 as pm
import theano.tensor as tt
from pymc3.math import logsumexp
from matplotlib import pyplot as plt
import sys
pm.distributions.continuous.HalfCauchy

# def mixture_density(w, mu, sd, x):
#     logp = tt.log(w) + pm.Normal.dist(mu, sd).logp(x)
#     return tt.sum(tt.exp(logp), axis=1)
#
# x = np.linspace(-5, 10, 100)[:, None]
# mu0 = np.array([-1., 2.6])
# sd0 = np.array([.5, 1.4])
# w0 =np.array([10, 60])
# yhat = mixture_density(w0, mu0, sd0, x).eval()
# y = yhat + np.random.randn(100)
# plt.plot(x, yhat);
# plt.scatter(x, y);
# plt.show()

# import numpy as np
# import pymc3 as pm
# import theano.tensor as tt
# import matplotlib.pyplot as plt
# from pymc3.math import logsumexp
#
# gaussCont = 1/np.sqrt(2*np.pi)
#
#
# def mixture_density_single(w, mu, sd, x):
#     logp = tt.log(w) + pm.Normal.dist(mu, sd).logp(x)
#     return tt.exp(logp)
#
#
# def mixture_density_mult(w, mu, sd, x):
#     #logp = tt.log(w) + pm.Normal.dist(mu, sd).logp(x)
#     return pm.Normal.dist(mu, sd).logp(x)#tt.sum(tt.exp(logp), axis=1)
#
#
# def mixture_density(w, mu, sd, x):
#     logp = tt.log(w) + pm.Normal.dist(mu, sd).logp(x)
#     return tt.sum(tt.exp(logp), axis=1)
#
#
# def gaussFunc(ind_params, a, mu, sigma):
#     x, z = ind_params
#     return a * gaussCont/sigma * np.exp(-((x - mu) * (x - mu)) / (2 * (sigma * sigma))) + z
#
#
# x = np.linspace(-5, 10, 100)[:, None]
# x_mine = np.linspace(-5, 10, 100)
# mu0 = np.array([-1., 2.6])
# sd0 = np.array([.5, 1.4])
# w0 = np.array([10, 60])
# y_orig = mixture_density(w0, mu0, sd0, x).eval()
# y_single = mixture_density_single(w0[0], mu0[0], sd0[0], x).eval()
#
# x_2d = np.atleast_2d(x_mine)
# y_mult = mixture_density_mult(w0[:, None], mu0[:, None], sd0[:, None], x_2d).eval()
#
# print('Aqui paro')
#
# # Plot data
# for i in np.arange(w0.size):
#     y_i = gaussFunc((x_mine, np.zeros(x.size)), w0[i], mu0[i], sd0[i])
#     plt.plot(x, y_i, label=f'Gaussian {i}')
#
# plt.plot(x, y_orig, label='Original function')
# plt.plot(x, y_mult, label='My function')
# plt.legend()
# plt.show()



# 2.28, 4.08
# 8.53, 17.05

# yhat = mixture_density(w0, mu0, sd0, x).eval()
# yhat2 = mixture_density([1,1], mu0, sd0, x).eval()
#
#
#
# plt.plot(x, yhat2)
#

# with pm.Model():
#     w = pm.HalfNormal('w', 10., shape=2)
#     mu = pm.Normal('mu', 0., 100., shape=2)
#     sd = pm.HalfCauchy('sd', 5., shape=2)
#     #     noise = pm.HalfCauchy('eps', 5.)
#     pm.Normal('obs', mixture_density(w, mu, sd, x), 1., observed=y)
#     trace = pm.sample(chains=2, cores=1)
#
# # pm.traceplot(trace)
# pm.plot_posterior(trace)


# import matplotlib.pyplot as plt
# import numpy as np
# import scipy.stats as st
# plt.style.use('seaborn-darkgrid')
# x = np.linspace(0, 3, 100)
# mus = [0., 0., 0.]
# sigmas = [.25, .5, 1.]
# for mu, sigma in zip(mus, sigmas):
#     pdf = st.lognorm.pdf(x, sigma, scale=np.exp(mu))
#     plt.plot(x, pdf, label=r'$\mu$ = {}, $\sigma$ = {}'.format(mu, sigma))
# plt.xlabel('x', fontsize=12)
# plt.ylabel('f(x)', fontsize=12)
# plt.legend(loc=1)
# plt.show()

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
plt.style.use('seaborn-darkgrid')
x = np.linspace(0, 5, 200)
for b in [0.5, 2.0, 5.0]:
    pdf = st.cauchy.pdf(x, scale=b)
    plt.plot(x, pdf, label=r'$\beta$ = {}'.format(b))
plt.xlabel('x', fontsize=12)
plt.ylabel('f(x)', fontsize=12)
plt.legend(loc=1)
plt.show()

# np.random.standard_cauchy()
#
# import matplotlib.pyplot as plt
# import numpy as np
# import scipy.stats as st
#
# plt.style.use('seaborn-darkgrid')
# x = np.linspace(-5000, 20000, 10000)
# alphas = [0., 0., 0., -2.]
# betas = [3, 1., 2., 1.]
# param0, param1 = 1000, 15000
# for a, b in zip(alphas, betas):
#     pdf = st.cauchy.pdf(x, loc=a, scale=b)
#     plt.plot(x, pdf * param0 + param1, label=r'$\alpha$ = {}, $\beta$ = {}'.format(a, b))
# plt.xlabel('x', fontsize=12)
# plt.ylabel('f(x)', fontsize=12)
# plt.legend(loc=1)
# plt.show()


from scipy.stats import cauchy
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(1, 1)



# Author: Jake VanderPlas
# License: BSD
#   The figure produced by this code is published in the textbook
#   "Statistics, Data Mining, and Machine Learning in Astronomy" (2013)
#   For more information, see http://astroML.github.com
#   To report a bug or issue, use the following forum:
#    https://groups.google.com/forum/#!forum/astroml-general
# import numpy as np
# from scipy.stats import cauchy
# from matplotlib import pyplot as plt
#
#
# #------------------------------------------------------------
# # Define the distribution parameters to be plotted
# gamma_values = [500, 1000, 2000]
# linestyles = ['-', '--', ':']
# mu = 0
# x = np.linspace(-5000, 20000, 10000)
# param0, param1 = 1000, 15000
# #------------------------------------------------------------
# # plot the distributions
# fig, ax = plt.subplots(figsize=(5, 3.75))
#
# for gamma, ls in zip(gamma_values, linestyles):
#     dist = cauchy(12000, gamma)
#     #y = 1 / (np.pi * gamma * (1 + np.power((x - mu)/gamma, 2)))
#     # plt.plot(x + param1, y * param0 , ls=ls, label=f'Mine {mu} {gamma}')
#     hist = 12000 + gamma * np.random.standard_cauchy(size=1000000)
#     hist = hist[(hist > 8000) & (hist < 20000)]
#     plt.hist(hist, bins=1000, normed=True)
#     plt.plot(x, dist.pdf(x), ls=ls, color='black', label=r'$\mu=%i,\ \gamma=%.1f$' % (mu, gamma))
#
# plt.xlabel('$x$')
# plt.ylabel(r'$p(x|\mu,\gamma)$')
# plt.title('Cauchy Distribution')
#
# plt.legend()
# plt.show()

