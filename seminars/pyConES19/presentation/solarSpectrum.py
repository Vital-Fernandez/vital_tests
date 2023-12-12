import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

background = np.array((43, 43, 43))/255.0
foreground = np.array((179, 199, 216))/255.0
red = np.array((43, 43, 43))/255.0
yellow = np.array((191, 144, 0))/255.0

data_folder = 'E:/Dropbox/Astrophysics/Seminars/PyConES_2019/'
file_address = 'astmg173.xls'
columns_names = ['Wavelength', 'Flux', 'Global_Tilt', 'Direct_circumsolar']

spec_df = pd.read_excel(data_folder+file_address, names=columns_names, delim_whitespace=True, skiprows=2)

idx_wave = spec_df['Wavelength'] < 3000
wavelength_file = spec_df['Wavelength'].values[idx_wave]

h = 6.626e-34
c = 3e8
T = 6000
k = 1.38066e-23
wavelength = np.arange(0.2, 3200) * 1e-9

p = 2*h*c*c/(wavelength**5)
b6000 = p/(np.exp(h*c/(wavelength*k*T))-1)
b6000 = (1e-9)*b6000
b6000 = b6000*(2.177e-5);
b6000 = b6000*np.pi
wavelength_nm = wavelength * 1e6

figConf = {'text.color': foreground,
            'figure.figsize': (16,8),
            'figure.facecolor':background,
            'axes.facecolor':background,
            'axes.edgecolor':foreground,
            'axes.labelcolor':foreground,
            'axes.labelsize':20,
            'xtick.color':foreground,
            'ytick.color':foreground,
            'xtick.labelsize':16,
            'ytick.labelsize':16,
            'legend.edgecolor':'inherit',
            'legend.facecolor':'inherit',
            'legend.fontsize':16,
            'legend.loc': "center right"}

matplotlib.rcParams.update(figConf)

scale = 1000
wave_blackBody = wavelength_nm * scale
wave_spectrum = spec_df['Wavelength'].values/1e3 * scale

fig, ax = plt.subplots()
ax.plot(wave_blackBody, b6000, label='Black body ~6000K', color=foreground)
ax.plot(wave_spectrum, spec_df['Flux'].values, label='Solar spectrum', color=yellow)

ax.set_xlabel('Photon Wavelength (nm)')
ax.set_ylabel('Number of photons')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# ax.tick_params(which='both', axis='x', bottom=False, top=False, labelbottom=False)
ax.tick_params(which='both', axis='y', left=False, labelleft=False)
ax.set_xlim(xmax=2500)
ax.legend()

# plt.show()
plt.savefig(data_folder+'sunSpectrum.png', tight=True, facecolor=background)







#
# fig, ax = plt.subplots()
# idx_wave = spec_df['Wavelength'] < 2000
# # ax.plot(spec_df['Wavelength'].values[idx_wave], spec_df['Flux'].values[idx_wave], label='Black body model at 2.725 K', color=foreground)
# ax.plot(wavelength, b6000)
#
# ax.set_xlabel('Photon Wavelength')
# ax.set_ylabel('Number of photons')
# ax.legend()
#
# plt.show()