from os import environ
environ["MKL_THREADING_LAYER"] = "GNU"

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import norm

# sns.set_style('white')
# sns.set_context('talk')

background = np.array((43, 43, 43))/255.0
foreground = np.array((179, 199, 216))/255.0
red = np.array((43, 43, 43))/255.0
yellow = np.array((191, 144, 0))/255.0

data_folder = 'E:/Dropbox/Astrophysics/Seminars/PyConES_2019/'

figConf = {'text.color': foreground,
            # 'figure.figsize': (16,10),
            'figure.facecolor':background,
            'axes.facecolor':background,
            'axes.edgecolor':foreground,
            'axes.labelcolor':foreground,
            # 'axes.labelsize':18,
            # 'xtick.labelsize':16,
            # 'ytick.labelsize':16,
            'xtick.color':foreground,
            'ytick.color':foreground,
            'legend.edgecolor':'inherit',
            'legend.facecolor':'inherit',
            # 'legend.fontsize':16,
             'legend.loc':"center right"}

matplotlib.rcParams.update(figConf)








np.random.seed(144) #122 124

data = np.random.randn(20)

def calc_posterior_analytical(data, x, mu_0, sigma_0):
    sigma = 1.
    n = len(data)
    mu_post = (mu_0 / sigma_0**2 + data.sum() / sigma**2) / (1. / sigma_0**2 + n / sigma**2)
    sigma_post = (1. / sigma_0**2 + n / sigma**2)**-1
    return norm(mu_post, np.sqrt(sigma_post)).pdf(x)

def sampler(data, samples=4, mu_init=.5, proposal_width=.5, plot=False, mu_prior_mu=0, mu_prior_sd=1.):
    mu_current = mu_init
    posterior = [mu_current]
    for i in range(samples):
        # suggest new position
        mu_proposal = norm(mu_current, proposal_width).rvs()

        # Compute likelihood by multiplying probabilities of each data point
        likelihood_current = norm(mu_current, 1).pdf(data).prod()
        likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()

        # Compute prior probability of current and proposed mu
        prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
        prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)

        p_current = likelihood_current * prior_current
        p_proposal = likelihood_proposal * prior_proposal

        # Accept proposal?
        p_accept = p_proposal / p_current

        # Usually would include prior probability, which we neglect here for simplicity
        random_factor = np.random.rand()
        accept = random_factor < p_accept

        if plot:
            plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accept, posterior, i, random_factor, p_accept)

        if accept:
            # Update position
            mu_current = mu_proposal

        posterior.append(mu_current)

    return posterior

# Function to display
def plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accepted, trace, i, random_factor, p_accept):

    from copy import copy
    trace = copy(trace)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, figsize=(16, 4))
    fig.suptitle('Iteration %i' % (i + 1))
    x = np.linspace(-3, 3, 5000)
    color = 'g' if accepted else 'r'

    # Plot prior
    prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
    prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)
    prior = norm(mu_prior_mu, mu_prior_sd).pdf(x)
    ax1.plot(x, prior)
    ax1.plot([mu_current] * 2, [0, prior_current], marker='o', color='b')
    ax1.plot([mu_proposal] * 2, [0, prior_proposal], marker='o', color=color)
    ax1.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    ax1.set(ylabel='Probability Density', title='current: prior(mu=%.2f) = %.2f\nproposal: prior(mu=%.2f) = %.2f' % (mu_current, prior_current, mu_proposal, prior_proposal))

    # Likelihood
    likelihood_current = norm(mu_current, 1).pdf(data).prod()
    likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()
    y = norm(loc=mu_proposal, scale=1).pdf(x)
    #y = np.array([norm(loc=i, scale=1).pdf(data).prod() for i in x])
    sns.distplot(data, kde=False, norm_hist=True, ax=ax2)
    ax2.plot(x, y, color=color)
    ax2.axvline(mu_current, color='b', linestyle='--', label='mu_current')
    ax2.axvline(mu_proposal, color=color, linestyle='--', label='mu_proposal')
    #ax2.title('Proposal {}'.format('accepted' if accepted else 'rejected'))
    ax2.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    ax2.set(title='likelihood(mu=%.2f) = %.2f\nlikelihood(mu=%.2f) = %.2f' % (mu_current, 1e14*likelihood_current, mu_proposal, 1e14*likelihood_proposal))

    # Posterior
    posterior_analytical = calc_posterior_analytical(data, x, mu_prior_mu, mu_prior_sd)
    ax3.plot(x, posterior_analytical)
    posterior_current = calc_posterior_analytical(data, mu_current, mu_prior_mu, mu_prior_sd)
    posterior_proposal = calc_posterior_analytical(data, mu_proposal, mu_prior_mu, mu_prior_sd)
    ax3.plot([mu_current] * 2, [0, posterior_current], marker='o', color='b')
    ax3.plot([mu_proposal] * 2, [0, posterior_proposal], marker='o', color=color)
    ax3.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    #x3.set(title=r'prior x likelihood $\propto$ posterior')
    ax3.set(title='posterior(mu=%.2f) = %.5f\nposterior(mu=%.2f) = %.5f' % (mu_current, posterior_current, mu_proposal, posterior_proposal))

    if accepted:
        trace.append(mu_proposal)
    else:
        trace.append(mu_current)
    ax4.plot(trace)
    operation = '>' if p_accept > random_factor else '<'
    title_label = 'Trace\nP(accept)={:.2E} {} P(random)={:.2E}'.format(p_accept, operation, random_factor)
    ax4.set(xlabel='iteration', ylabel='mu', title=title_label)
    plt.tight_layout()
    # plt.legend()
    # plt.show()
    plt.savefig(data_folder + f'iteration_{i}.png', bbox_inches='tight', facecolor=background)


# Plot of iteration steps
sampler(data, samples=16, mu_init=-1., plot=True)

# # # Plot of good trace
# posterior = sampler(data, samples=5000, mu_init=1.)
# fig, ax = plt.subplots()
# ax.plot(posterior)
# _ = ax.set(xlabel='sample', ylabel='$T_{\mu}$');
# # plt.show()
#
# # Plot too narow proposal width
# posterior_small = sampler(data, samples=5000, mu_init=1., proposal_width=.01)
# fig, ax = plt.subplots()
# ax.plot(posterior_small);
# _ = ax.set(xlabel='sample', ylabel='$T_{\mu}$');
# # plt.show()
#
# # Plot too wide proposal width proposal width
# posterior_large = sampler(data, samples=5000, mu_init=1., proposal_width=3.)
# fig, ax = plt.subplots()
# ax.plot(posterior_large); plt.xlabel('sample'); plt.ylabel('mu');
# _ = ax.set(xlabel='sample', ylabel='$T_{\mu}$');
# # plt.show()
#
#
# from pymc3.stats import autocorr
# lags = np.arange(1, 100)
# fig, ax = plt.subplots()
# ax.plot(lags, [autocorr(posterior_large, l) for l in lags], label='Wide proposal width')
# ax.plot(lags, [autocorr(posterior_small, l) for l in lags], label='Narrow proposal width')
# ax.plot(lags, [autocorr(posterior, l) for l in lags], label='Medium proposal width')
# ax.legend(loc=0)
# _ = ax.set(xlabel='lag', ylabel='autocorrelation', ylim=(-.1, 1))
#
# plt.show()