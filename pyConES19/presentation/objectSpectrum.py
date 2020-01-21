import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.io import fits

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

figConf = {'text.color': foreground,
            'figure.figsize': (16,10),
            'figure.facecolor':background,
            'axes.facecolor':background,
            'axes.edgecolor':foreground,
            'axes.labelcolor':foreground,
            'axes.labelsize':18,
            'xtick.labelsize':16,
            'ytick.labelsize':16,
            'xtick.color':foreground,
            'ytick.color':foreground,
            'legend.edgecolor':'inherit',
            'legend.facecolor':'inherit',
            'legend.fontsize':16,
             'legend.loc':"center right"}

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

plt.savefig(data_folder+'objectSpectrum.png', bbox_inches='tight', facecolor=background)


# plt.show()
