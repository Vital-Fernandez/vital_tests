import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


background = np.array((43, 43, 43))/255.0
foreground = np.array((179, 199, 216))/255.0
red = np.array((43, 43, 43))/255.0
yellow = np.array((191, 144, 0))/255.0

data_folder = 'E:/Dropbox/Astrophysics/Seminars/PyConES_2019/'
spec_folder = 'E:/Dropbox/Astrophysics/Data/WHT_observations/objects/8/'
spec_Name = '8_WHT.fits'

data_array, Header_0 = fits.getdata(spec_folder + spec_Name, header=True)

wavelength = data_array['Wave'] / 10
flux = data_array['Int'] * 1e14

idx_region = (432 < wavelength) & (wavelength < 439)

figConf = {'text.color': foreground,
            'figure.figsize': (16,5),
            'figure.facecolor': background,
            'axes.facecolor': background,
            'axes.edgecolor': foreground,
            'axes.labelcolor': foreground,
            'axes.labelsize': 28,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'xtick.color': foreground,
            'ytick.color': foreground,
            'legend.edgecolor': 'inherit',
            'legend.facecolor': 'inherit',
            'legend.fontsize': 16,
             'legend.loc': "center right"}

matplotlib.rcParams.update(figConf)

fig, ax = plt.subplots()
fig.set_facecolor(background)
ax.plot(wavelength, flux, label='', color=foreground)
ax.set_xlabel('Photon Wavelength (nm)')
ax.set_ylabel('Number of photons')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# ax.tick_params(which='both', axis='x', bottom=False, top=False, labelbottom=False)
ax.tick_params(which='both', axis='y', left=False, labelleft=False)
# ax.legend()



# Plot zoomed region
axins = zoomed_inset_axes(ax, zoom=6, loc=6)
axins.step(wavelength[idx_region], flux[idx_region], where='mid', color=foreground)
axins.tick_params(which='both', axis='y', left=False, labelleft=False)
axins.tick_params(which='both', axis='x', bottom=False, top=False, labelbottom=False)
# axins.set_xlim(7550, 7700)
# axins.set_ylim(-0.2e-16, 2.0e-16)
mark_inset(ax, axins, loc1=3, loc2=4, fc="none", ec="0.5")

plt.savefig(data_folder+'detailObjectSpectrum.png', bbox_inches='tight', facecolor=background)
# plt.show()
