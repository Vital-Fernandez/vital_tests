# import numpy as np
# from pymc3 import *
# import matplotlib.pyplot as plt
#
# # set random seed for reproducibility
# np.random.seed(12345)
#
# x = np.arange(5,400,10)*1e3
#
# # Parameters for gaussian
# amp_true = 0.2
# size_true = 1.8
# ps_true = 0.1
#
# #Gaussian function
# gauss = lambda x,amp,size,ps: amp*np.exp(-1*(np.pi**2/(3600.*180.)*size*x)**2/(4.*np.log(2.)))+ps
# f_true = gauss(x=x,amp=amp_true, size=size_true, ps=ps_true )
#
# # add noise to the data points
# noise = np.random.normal(size=len(x)) * .02
# f = f_true + noise
# f_error = np.ones_like(f_true)*0.05*f.max()
#
# with Model() as model3:
#     amp = Uniform('amp', 0.05, 0.4, testval= 0.15)
#     size = Uniform('size', 0.5, 2.5, testval= 1.0)
#     ps = Normal('ps', 0.13, 40, testval=0.15)
#
#     gauss = Deterministic('gauss', amp*np.exp(-1*(np.pi**2*size*x/(3600.*180.))**2/(4.*np.log(2.)))+ps)
#
#     y =Normal('y', mu=gauss, tau=1.0/f_error**2, observed=f)
#
#     start=find_MAP()
#     step=NUTS()
#     trace=sample(2000, start=start)
#
# # extract and plot results
# y_min = np.percentile(trace.gauss,2.5,axis=0)
# y_max = np.percentile(trace.gauss,97.5,axis=0)
# y_fit = np.percentile(trace.gauss,50,axis=0)
# plt.plot(x, f_true, 'b', marker='None', ls='-', lw=1, label='True')
# plt.errorbar(x, f, yerr=f_error, color='r', marker='.', ls='None', label='Observed')
# plt.plot(x, y_fit, 'k', marker='+', ls='None', ms=5, mew=1, label='Fit')
# plt.fill_between(x, y_min, y_max, color='0.5', alpha=0.5)
# plt.legend()
# plt.show()

from matplotlib import pyplot as plt
import numpy as np
import pymc3 as pm
import scipy as sp
from theano import tensor as T

np.random.seed(462233) # from random.org

N = 20
K = 30

alpha = 2.
P0 = sp.stats.norm

beta = sp.stats.beta.rvs(1, alpha, size=(N, K))
w = np.empty_like(beta)
w[:, 0] = beta[:, 0]
w[:, 1:] = beta[:, 1:] * (1 - beta[:, :-1]).cumprod(axis=1)

omega = P0.rvs(size=(N, K))

x_plot = np.linspace(-3, 3, 200)

sample_cdfs = (w[..., np.newaxis] * np.less.outer(omega, x_plot)).sum(axis=1)


# # ---------------------------------------------------------
# fig, ax = plt.subplots(figsize=(8, 6))
#
# ax.plot(x_plot, sample_cdfs[0], c='gray', alpha=0.75,
#         label='DP sample CDFs');
# ax.plot(x_plot, sample_cdfs[1:].T, c='gray', alpha=0.75)
# ax.plot(x_plot, P0.cdf(x_plot), c='k', label='Base CDF')
#
# ax.set_title(r'$\alpha = {}$'.format(alpha))
# ax.legend(loc=2)
# plt.show()
#
# #----------------------------------------------------------
#
# fig, (l_ax, r_ax) = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(16, 6))
#
# K = 50
# alpha = 10.
#
# beta = sp.stats.beta.rvs(1, alpha, size=(N, K))
# w = np.empty_like(beta)
# w[:, 0] = beta[:, 0]
# w[:, 1:] = beta[:, 1:] * (1 - beta[:, :-1]).cumprod(axis=1)
#
# omega = P0.rvs(size=(N, K))
#
# sample_cdfs = (w[..., np.newaxis] * np.less.outer(omega, x_plot)).sum(axis=1)
#
# l_ax.plot(x_plot, sample_cdfs[0], c='gray', alpha=0.75,
#           label='DP sample CDFs');
# l_ax.plot(x_plot, sample_cdfs[1:].T, c='gray', alpha=0.75);
# l_ax.plot(x_plot, P0.cdf(x_plot), c='k', label='Base CDF');
#
# l_ax.set_title(r'$\alpha = {}$'.format(alpha));
# l_ax.legend(loc=2);
#
# K = 200
# alpha = 50.
#
# beta = sp.stats.beta.rvs(1, alpha, size=(N, K))
# w = np.empty_like(beta)
# w[:, 0] = beta[:, 0]
# w[:, 1:] = beta[:, 1:] * (1 - beta[:, :-1]).cumprod(axis=1)
#
# omega = P0.rvs(size=(N, K))
#
# sample_cdfs = (w[..., np.newaxis] * np.less.outer(omega, x_plot)).sum(axis=1)
#
# r_ax.plot(x_plot, sample_cdfs[0], c='gray', alpha=0.75,
#           label='DP sample CDFs');
# r_ax.plot(x_plot, sample_cdfs[1:].T, c='gray', alpha=0.75);
# r_ax.plot(x_plot, P0.cdf(x_plot), c='k', label='Base CDF');
#
# r_ax.set_title(r'$\alpha = {}$'.format(alpha));
# r_ax.legend(loc=2);
# plt.show()
# #----------------------------------------------------------

N = 5
K = 30

alpha = 2
P0 = sp.stats.norm
f = lambda x, theta: sp.stats.norm.pdf(x, theta, 0.3)

beta = sp.stats.beta.rvs(1, alpha, size=(N, K))
w = np.empty_like(beta)
w[:, 0] = beta[:, 0]
w[:, 1:] = beta[:, 1:] * (1 - beta[:, :-1]).cumprod(axis=1)

theta = P0.rvs(size=(N, K))

dpm_pdf_components = f(x_plot[np.newaxis, np.newaxis, :], theta[..., np.newaxis])
dpm_pdfs = (w[..., np.newaxis] * dpm_pdf_components).sum(axis=1)

fig, ax = plt.subplots(figsize=(8, 6))

ax.plot(x_plot, dpm_pdfs.T, c='gray');

ax.set_yticklabels([]);

plt.show()

fig, ax = plt.subplots(figsize=(8, 6))

ix = 1

ax.plot(x_plot, dpm_pdfs[ix], c='k', label='Density');
ax.plot(x_plot, (w[..., np.newaxis] * dpm_pdf_components)[ix, 0],
        '--', c='k', label='Mixture components (weighted)');
ax.plot(x_plot, (w[..., np.newaxis] * dpm_pdf_components)[ix].T,
        '--', c='k');

ax.set_yticklabels([]);
ax.legend(loc=1);

plt.show()