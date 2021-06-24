import numpy as np
import pyneb as pn
from matplotlib import pyplot as plt, rcParams
from astro.data.muse.common_methods import STANDARD_AXES, DARK_PLOT, background_color, foreground_color

# %matplotlib inline

defaultConf = DARK_PLOT.copy()
rcParams.update(defaultConf)

size_dist = 10000
Te_dist = np.random.normal(10000, 5000, size_dist)
ne_dist = np.random.normal(250, 150, size_dist)

HI = pn.RecAtom('H', 1)

H_alphaBeta = HI.getEmissivity(tem=Te_dist, den=ne_dist, wave=6563, product=False)/HI.getEmissivity(tem=Te_dist, den=ne_dist, wave=4861, product=False)

print(H_alphaBeta.size)

fig = plt.figure()
ax = fig.add_subplot()
ax.hist(H_alphaBeta, bins=50, histtype='step')
ax.set_xlabel(r'$H\alpha/H\beta$ emissivity ratio')
ax.set_title(r'$H\alpha/H\beta$ emissivity ratio where $Te (\mu=10000K,\,\sigma=5000K$) and $ne(\mu=250cm^{-3},\,\sigma=150cm^{-3}$)')
plt.savefig('/home/vital/HalphaBeta_variation.png', dpi=300, bbox_inches='tight')
