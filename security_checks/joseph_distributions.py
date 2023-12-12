import numpy as np
import pandas as pd
from pathlib import Path
import lime
from matplotlib import pyplot as plt, rc_context


PLOT_FORMAT = {'figure.figsize': (8, 4), 'axes.titlesize': 18, 'axes.labelsize': 25, 'legend.fontsize': 10,
               'xtick.labelsize': 10, 'ytick.labelsize': 10, 'font.family': 'Times New Roman', 'mathtext.fontset':'cm',
               "figure.dpi": 200}

file_path = Path(f'D:/Dropbox/Aplicaciones/joseph_data/edited data_v3.xlsx')

df = lime.load_log(file_path, page='Sheet1')

iq_mean, iq_std = 100, 15
iq_dist = np.random.normal(100, 15, size=df.index.size)

with rc_context(PLOT_FORMAT):

    fig, ax = plt.subplots()

    array1 = df['1st IQ score'].to_numpy()
    array2 = df['2nd IQ scores'].to_numpy()

    ax.hist(array1, histtype='stepfilled', alpha=.5, density=False, label=r'$1^{st}$ IQ score')
    ax.hist(array2, histtype='stepfilled', alpha=.5, density=False, label=r'$2^{nd}$ IQ score')
    # ax.hist(iq_dist, histtype='step', alpha=.5, density=False, color='black', label=r'Literature $\mathcal{N}(100, 15)$')

    ax.update({'xlabel': 'Intelligence Quotient', 'ylabel': 'Count (111 total)',
               'title': 'IQ distribution from measurements'})

    ax.axvline(iq_mean, color='black')
    ax.axvspan(iq_mean - iq_std, iq_mean + iq_std, alpha=0.5, color='grey')

    # ax.hist(trace, histtype='stepfilled', bins=35, alpha=.7, color=cmap(colorNorm(colorDict[ion])), density=False)
    ax.legend()

    plt.show()
