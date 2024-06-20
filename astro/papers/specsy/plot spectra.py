import lime
import numpy as np
from matplotlib import pyplot as plt, rc_context
from pathlib import Path
from lime.tools import au
lime.theme.set_style('dark')
fig_cfg = lime.theme.fig_defaults()
ax_cfg = lime.theme.ax_defaults(units_flux=au.Unit('FLAM'), units_wave=au.Unit('AA'), norm_flux=1, user_ax=None)
ax_cfg['ylabel'] = 'Normalized flux'

fig_cfg['figure.figsize'] = (5,8)
fig_cfg['legend.fontsize'] = 15
fig_cfg['axes.labelsize'] = 20

fig_cfg['ytick.left'] = False
fig_cfg['ytick.labelleft'] = False
fig_cfg['ytick.right'] = False
fig_cfg['ytick.labelright'] = False

data_folder = Path('/home/vital/Astrodata/Z004')
spectra_list = ['/home/vital/Astrodata/Z004/SSP-CHA-stellar_Z0.004_logt5.00.dat',
                '/home/vital/Astrodata/Z004/SSP-CHA-stellar_Z0.004_logt6.00.dat',
                '/home/vital/Astrodata/Z004/SSP-CHA-stellar_Z0.004_logt7.00.dat',
                '/home/vital/Astrodata/Z004/SSP-CHA-stellar_Z0.004_logt8.00.dat',
                '/home/vital/Astrodata/Z004/SSP-CHA-stellar_Z0.004_logt9.00.dat']

w0, wf = 3600, 7000
wN1, wN2 = 6000, 6200
label_list = ['log(age) = 5', 'log(age) = 6', 'log(age) = 7', 'log(age) = 8', 'log(age) = 9']

with rc_context(fig_cfg):

    fig, ax = plt.subplots()

    for i, address in enumerate(spectra_list):

        wave, flux = np.loadtxt(address, unpack=True, skiprows=1)
        idcs = np.searchsorted(wave, (w0, wf))

        idcs_norm = np.searchsorted(wave, (wN1, wN2))
        flux_norm = np.sum(flux[idcs_norm[0]:idcs_norm[1]])

        wave, flux = wave[idcs[0]:idcs[1]], flux[idcs[0]:idcs[1]]/flux_norm
        ax.step(wave, flux, linewidth=0.5, label=label_list[i])

    ax.update(ax_cfg)
    ax.legend()
    plt.tight_layout()
    #plt.show()
    plt.savefig('/home/vital/Desktop/SSP_popstar.png')
