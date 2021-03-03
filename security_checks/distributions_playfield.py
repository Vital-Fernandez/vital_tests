import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import truncnorm, norm
from astropy import modeling
import scipy.stats as st


params_prior_conf = {'Teff_prior': ['Uniform', 30000.0, 90000.0, 1.0, 0.0],
                     'logU_prior': ['Uniform', -4.0, -1.5, 1.0, 0.0],
                     'T_low_prior': ['Normal', 15000.0, 5000.0, 1.0, 0.0],
                     'T_high_prior': ['Normal', 15000.0, 5000.0, 1.0, 0.0],
                     'n_e_prior': ['HalfCauchy', 2.0, 0.0, 200.0, 0.0],
                     'cHbeta_prior': ['HalfCauchy', 2.0, 0.0, 1.0, 0.0],
                     'tau_prior': ['Lognormal', 0, 1.0, 1.0, 0.0],
                      'He1_prior': ['Normal', 0.0, 3.0, 1.0, 0.0], 'He2_prior': ['Normal', 0.0, 3.0, 1.0, 0.0],
                     'Ar3_prior': ['Normal', 5.0, 5.0, 1.0, 0.0], 'Ar4_prior': ['Normal', 5.0, 5.0, 1.0, 0.0],
                     'N2_prior': ['Normal', 5.0, 5.0, 1.0, 0.0], 'O2_prior': ['Normal', 5.0, 5.0, 1.0, 0.0],
                     'O3_prior': ['Normal', 5.0, 5.0, 1.0, 0.0], 'S2_prior': ['Normal', 5.0, 5.0, 1.0, 0.0],
                     'S3_prior': ['Normal', 5.0, 5.0, 1.0, 0.0], 'Cl3_prior': ['Normal', 5.0, 5.0, 1.0, 0.0],
                     'Ne3_prior': ['Normal', 5.0, 5.0, 1.0, 0.0], 'Fe3_prior': ['Normal', 5.0, 5.0, 1.0, 0.0]}

# def gauss_curve(x, amp_g, mu_g, sigma_g):
#     y = 1 / (sigma_g * np.sqrt(2 * np.pi)) * np.exp(- (x - mu_g) ** 2 / (2 * sigma_g ** 2))
#     return y
#
# def trunc_limits(lower_limit, upper_limit, mu_g, sigma_g):
#     return (lower_limit - mu_g) / sigma_g, (upper_limit - mu_g) / sigma_g
#
# n_steps = 100000
#
# mu, sigma = 1, 0.8 # mean and standard deviation
# s = np.random.normal(mu, sigma, n_steps)
#
# a_trunc, b_trunc = trunc_limits(0.0, np.infty, mu, sigma)
# s_trunc = truncnorm.rvs(a_trunc, b_trunc, loc=mu, scale=sigma, size=n_steps)
#
# print(abs(mu - np.mean(s)))
# print(abs(sigma - np.std(s, ddof=1)))
# print(a_trunc, b_trunc)
# count, bins, ignored = plt.hist(s, 30, density=True, range=(-0.5,5.0))
# count_trunc, bins_trunc, ignored_trunc = plt.hist(s_trunc, 30, density=True, alpha=0.5, range=(-0.5,5.0))
#
# fitter = modeling.fitting.LevMarLSQFitter()
# model = modeling.models.Gaussian1D()   # depending on the data you need to give some initial values
# fitted_model = fitter(model, bins, gauss_curve(x=bins, amp_g=0, mu_g=mu, sigma_g=sigma))
# fitted_model_trunc = fitter(model, bins_trunc, truncnorm.pdf(bins_trunc, a_trunc, b_trunc, loc=mu, scale=sigma))
#
# plt.plot(bins, gauss_curve(x=bins, amp_g=0, mu_g=mu, sigma_g=sigma), linewidth=2)
# plt.plot(bins_trunc, gauss_curve(x=bins_trunc, amp_g=0, mu_g=mu, sigma_g=sigma), label='truncnorm pdf')
# plt.show()

# import numpy as np
# from scipy.stats import truncnorm
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(1, 1)
#
# a, b = 0.1, 2
# mean, var, skew, kurt = truncnorm.stats(a, b, moments='mvsk')
# print(mean, var, skew, kurt)
#
# trunc_limits()
#
# x = np.linspace(truncnorm.ppf(0.01, a, b), truncnorm.ppf(0.99, a, b), 100)
# ax.plot(x, truncnorm.pdf(x, a, b), 'r-', lw=5, alpha=0.6, label='truncnorm pdf')
#
# rv = truncnorm(a, b)
# ax.plot(x, rv.pdf(x), 'k-', lw=2, label='frozen pdf')
# vals = truncnorm.ppf([0.001, 0.5, 0.999], a, b)
# np.allclose([0.001, 0.5, 0.999], truncnorm.cdf(vals, a, b))
#
# r = truncnorm.rvs(a, b, size=1000)
#
# ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
# ax.legend(loc='best', frameon=False)
# plt.show()



# x = np.linspace(-5000, 25000, 10000)
# alphas = [0., 0., 0., -2.]
# betas = [.5, 1., 2., 1.]
# for a, b in zip(alphas, betas):
#     pdf = st.cauchy.pdf(x, loc=a, scale=b)
#     plt.plot(x, pdf, label=r'$\alpha$ = {}, $\beta$ = {}'.format(a, b))
# plt.xlabel('x', fontsize=12)
# plt.ylabel('f(x)', fontsize=12)
# plt.legend(loc=1)
# plt.show()
