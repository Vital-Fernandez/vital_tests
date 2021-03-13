import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path

file_address = Path('D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/data/Papaderos_Full/6March2021/Base.BC03.All')

columns = ['spec-file', 'age', 'Z', 'code', 'Mstar', 'YAV', 'a/Fe']
basesDF = pd.read_csv(file_address, skiprows=1, header=None, delim_whitespace=True, names=columns, comment='#')

# Plot format
size_dict = {'figure.figsize': (18, 6), 'axes.titlesize': 16, 'axes.labelsize': 18, 'xtick.labelsize': 14,
             'ytick.labelsize': 14, 'legend.fontsize': 18}
rcParams.update(size_dict)

fig, ax = plt.subplots()
ax.scatter(np.log10(basesDF.age.values), basesDF.Z.values)
ax.update({'xlabel': 'log(Age) (yr)',
           'ylabel': r'Metallicity $(Z_\star)$',
           'title': f'{basesDF.index.size} SSPs Bruzual & Charlot 2003 library'})
plt.tight_layout()
plt.show()

# # Plot the results
# fig, ax = plt.subplots()
# ax.step(lm.wave, lm.flux, label='Object spectrum')
# ax.scatter(specWave, specFlux, label='Line', color='tab:orange')
# ax.plot(wave_resample, hmc_curve + cont_resample, label='HMC fitting',  color='tab:red')
# ax.legend()
# ax.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'Gaussian fitting'})
# plt.show()