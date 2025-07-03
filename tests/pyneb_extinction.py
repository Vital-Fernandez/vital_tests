import pyneb as pn
import numpy as np
from scipy.stats import truncnorm


def truncated_positive_normal(mean, std, n_steps=10000):
    """
    Generate samples from a normal distribution truncated to positive values.

    Parameters:
        mean (float): Mean of the original normal distribution.
        std (float): Standard deviation of the original normal distribution.
        size (int): Number of samples to generate (default 10,000).

    Returns:
        np.ndarray: Array of truncated normal samples > 0.
    """
    a, b = (0 - mean) / std, np.inf  # lower and upper bounds in standard normal units
    return truncnorm.rvs(a, b, loc=mean, scale=std, size=n_steps)


temp, den = 10000, 10000000
size = 10000

H1 = pn.RecAtom('H', 1)
O3 = pn.Atom('O', 3)

# arrays
Hbeta_arr = truncated_positive_normal(172.0, 24, size)
Hgamma_arr = truncated_positive_normal(55.6, 17.2, size)

emis_ratio = H1.getEmissivity(tem=temp, den=den, wave=4340)/H1.getEmissivity(tem=temp, den=den, wave=4861)
obs_ratio = Hgamma_arr/Hbeta_arr
print('emissivity ratio', emis_ratio)
print('observed ratio', obs_ratio)

rc = pn.RedCorr(R_V=3.4, law='G03 LMC')
rc.setCorr(obs_over_theo=obs_ratio/emis_ratio, wave1=4340., wave2=4861.)
print(np.nanmean(rc.E_BV), np.nanstd(rc.E_BV))

cHbeta_arr = truncated_positive_normal(np.nanmean(rc.cHbeta), np.nanstd(rc.cHbeta), )
print(rc.getCorrHb(4363))
print(rc.getErrCorrHb(4363, err_E_BV=np.nanstd(rc.E_BV)))