# Author: Jake VanderPlas <vanderplas@astro.washington.edu>
# License: BSD
#   The figure is an example from astroML: see http://astroML.github.com
from matplotlib import pyplot as plt
from astroML.datasets import fetch_sdss_filter, fetch_vega_spectrum
import numpy as np
from pathlib import Path
import lime

output_folder = Path(f'D:/Dropbox/Astrophysics/Seminars/Umich_introduction_2023')

spec_manga_file = 'D:/Pycharm Projects/lime/tests/data_tests/manga_spaxel.txt'
wave_array, flux_array, err_array = np.loadtxt(spec_manga_file, unpack=True)

conf = {'figure.figsize': (10, 5),
          'axes.titlesize': 16,
          'axes.labelsize': 16,
          'legend.fontsize': 16,
          'xtick.labelsize': 14,
          'ytick.labelsize': 14}


spec = lime.Spectrum(wave_array, flux_array, err_array, 0.0475, norm_flux=1e-18)
spec.plot.spectrum(rest_frame=True, fig_cfg=conf, log_scale=True, output_address=output_folder/f'example_spectrum_log.png')

wave_rest = wave_array/(1 + 0.0475)
flux_array = flux_array/np.median(flux_array)/10
redshift = 0.4
#------------------------------------------------------------
# Set up figure and axes
fig = plt.figure(figsize=(15, 8))
ax = fig.add_subplot(111)

#----------------------------------------------------------------------
# Fetch and plot the Vega spectrum
spec = fetch_vega_spectrum()
lam = spec[0]
spectrum = spec[1] / 2.1 / spec[1].max()
# ax.plot(lam, spectrum, '-k', lw=2)
ax.plot(wave_rest * (1 + redshift), flux_array, '-k', lw=2)

#------------------------------------------------------------
# Fetch and plot the five filters
text_kwargs = dict(fontsize=20, ha='center', va='center', alpha=0.5)

for f, c, loc in zip('ugriz', 'bgrmk', [3500, 4600, 6100, 7500, 8800]):
    data = fetch_sdss_filter(f)
    ax.fill(data[0], data[1]*100, ec=c, fc=c, alpha=0.4)
    ax.text(loc, 5, f, color=c, **text_kwargs)

ax.set_xlim(3000, 11000)

ax.set_title(f'z = {redshift}',  fontsize = 40)
ax.set_xlabel('Wavelength (Angstroms)', fontsize = 20)

# plt.show()
plt.savefig(output_folder/f'sdss_redshift_{redshift}.png', bbox_inches='tight')
