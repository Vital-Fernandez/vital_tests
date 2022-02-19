import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from astro.papers.gtc_greenpeas.common_methods import red_corr_HalphaHbeta_ratio
from astro.data.muse.common_methods import STANDARD_PLOT, DARK_PLOT, background_color
from uncertainties import ufloat
from uncertainties import ufloat, unumpy
from src.specsiser.components.line_tools import STANDARD_PLOT, STANDARD_AXES
from matplotlib import pyplot as plt, rcParams

objList = ['GP030321', 'GP101157', 'GP121903', 'GP004054', 'GP113303', 'GP232539']

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, group_variables=False)
data_folder = Path(obsData['file_information']['data_folder'])

NO_ratio = np.array([-0.85, -1.40, -0.96, -1.06, -1.04, -0.99])
NO_ratio_err = np.array([0.06, 0.11, 0.04, 0.04, 0.08, 0.11])

OH_abund = np.array([7.76, 8.04, 7.94, 7.98, 7.91, 7.81])
OH_abund_err = np.array([0.02, 0.02, 0.02, 0.06, 0.1, 0.14])

OH_axis_log = np.linspace(7.50, 9.0, 50)

OH_axis = np.power(10, OH_axis_log - 12)
nicholls_relation = np.log10(np.power(10, -1.732) + np.power(10, np.log10(OH_axis) + 2.19))

x_gav_18, y_gav_18 = np.loadtxt(data_folder/'Gavilan_1_18kpc.csv', delimiter=',', unpack=True)
x_gav_8, y_gav_8 = np.loadtxt(data_folder/'Gavilan_2_8_kpc.csv', delimiter=',', unpack=True)
x_gav_4, y_gav_4 = np.loadtxt(data_folder/'Gavilan_3_4_kpc.csv', delimiter=',', unpack=True)
x_gav_2, y_gav_2 = np.loadtxt(data_folder/'Gavilan_4_2_kpc.csv', delimiter=',', unpack=True)

# Plot set up
defaultConf = STANDARD_PLOT.copy()

defaultConf['axes.labelsize'] = 20

rcParams.update(defaultConf)

fig, ax = plt.subplots(figsize=(12, 8))

ax.errorbar(OH_abund, NO_ratio, xerr=OH_abund_err, yerr=NO_ratio_err, fmt='o')
ax.plot(OH_axis_log, nicholls_relation, label='Nicholls et al. (2017)')
ax.plot(x_gav_2, y_gav_2, label='R = 2 Kpc (Gavilan et al. 2006)', linestyle=(0, (3, 1, 1, 1, 1, 1)))
ax.plot(x_gav_4, y_gav_4, label='R = 4 Kpc (Gavilan et al. 2006)', linestyle=':')
ax.plot(x_gav_8, y_gav_8, label='R = 8 Kpc (Gavilan et al. 2006)', linestyle='--')
ax.plot(x_gav_18, y_gav_18, label='R = 18 Kpc (Gavilan et al. 2006)', linestyle='dashdot')

ax.update({'xlabel': r'$12+log \left( \frac{O}{H} \right)$', 'ylabel': r'$log\left( \frac{N}{O}\right)$'})

ax.legend(loc='best')

ax.set_ylim(-2.0, 0.0)
ax.set_xlim(7.5, 8.5)

plt.show()

