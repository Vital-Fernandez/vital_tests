import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

data_folder = 'E:/Dropbox/Astrophysics/Seminars/PyConES_2019/'
background = np.array((43, 43, 43))/255.0
foreground = np.array((179, 199, 216))/255.0
file_name = 'firas_monopole_spec_v1.txt'

columns_names = ['frequency', 'flux', 'residual', 'uncertainty', 'ModeledSpec']

spec_df = pd.read_csv(data_folder + file_name, names=columns_names, delim_whitespace=True, skiprows=18)

figConf = {'text.color': foreground,
            'figure.figsize': (16,8),
            'figure.facecolor':background,
            'axes.facecolor':background,
            'axes.edgecolor':foreground,
            'axes.labelcolor':foreground,
            'axes.labelsize':18,
            'xtick.color':foreground,
            'ytick.color':foreground,
            'legend.edgecolor':'inherit',
            'legend.facecolor':'inherit',
            'legend.fontsize':16,
             'legend.loc':"center right"}

matplotlib.rcParams.update(figConf)

fig, ax = plt.subplots()
fig.set_facecolor(background)
ax.scatter(spec_df['frequency'], spec_df['flux'], label='Cosmic Microwave Background photons', color=foreground)
ax.plot(spec_df['frequency'], spec_df['flux'], label='Black body model at 2.725 K', color=foreground)
ax.set_xlabel('Photon Wavelength')
ax.set_ylabel('Number of photons')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(which='both', axis='x', bottom=False, top=False, labelbottom=False)
ax.tick_params(which='both', axis='y', left=False, labelleft=False)
ax.legend()

plt.savefig(data_folder+'bigBangSpectrum.png', tight=True, facecolor=background)
# plt.show()

# .tick_params(
#     axis='x',          # changes apply to the x-axis
#           # both major and minor ticks are affected
#           # ticks along the bottom edge are off
#              # ticks along the top edge are off
#     labelbottom=False) # labels along the bottom edge are off


# ax.spines['bottom'].set_color(foreground)
# ax.spines['top'].set_color(foreground)
# ax.xaxis.label.set_color(foreground)
# ax.yaxis.label.set_color(foreground)
# ax.tick_params(axis='x', colors=foreground)
# ax.tick_params(axis='y', colors=foreground)

