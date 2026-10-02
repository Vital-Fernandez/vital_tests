import sesamme.models as models
import sesamme.mcmc as stats
import sesamme.vis as vis
import emcee
import numpy as np
import corner
import lime
from astropy.table import Table
import astropy.units as u
### Data visualization
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm, ListedColormap, BoundaryNorm, Normalize
import matplotlib.patches as patches
from matplotlib.collections import LineCollection
import corner
from IPython.display import display, Math
from scipy.optimize import minimize
from specsy.models.stellar import StellarBinaries
from pathlib import Path
from astro.stsci.tools import get_masked_bands, CORNER_CFG


def log_likelihood_model(theta, wl, y, yerr, mask, model_cube, ion_table, add_nebular):
    y_model = models.get_model(theta, wl, model_cube, ion_table, add_nebular)

    masked_spec = np.array((y[mask] - y_model[mask]) / yerr[mask])
    masked_err = np.array(np.sqrt(2) * np.sqrt(np.pi) * yerr[mask])

    resid = -0.5 * (np.dot(masked_spec, masked_spec) + np.log(np.dot(masked_err, masked_err)))

    return resid


models.set_ext_law('Calzetti')

# ---get best initial parameters------#
# --we do this by minimizing the likelihood function that used by SESAMME

# ---set the distance in Mpc to your target to convert to luminosity
distance = 7.510  #

spec_pname = "/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs/NGC4861_pyStarburst99_spectrum_Lsun.txt"
spec = lime.Spectrum.from_file(spec_pname, instrument='text')
# spec.plot.spectrum()
wl = spec.wave.data
lum = spec.flux.data
lum_err = spec.err_flux.data
mask = spec.flux.mask

# read in data and models first:
cube_path = Path("/home/vital/Downloads/pySB99_SESAMME_Cube_v4/pySB99_SSP_Grid_v4.fits")
ssp_cube = StellarBinaries.from_fits(cube_path, source='pystarburst99')

obj_ssp_pname = f'./NGC4861_ssp_cube.fits'
ssp_cube.to_fits(obj_ssp_pname, disp_intvl=wl)
obj_cube = StellarBinaries.from_fits(obj_ssp_pname, source='pystarburst99')
modelcube = models.load_ssp_cube(obj_ssp_pname)

bin_wave, bin_flux = obj_cube.get_spectrum(6.046539389436416, 0.00040)
spec.plot.spectrum(in_fig=None)
spec.plot.ax.step(bin_wave, bin_flux, label='SSP log(age) = 6.046, z = 0.00040')
spec.plot.ax.legend()
spec.plot.show()


# modelcube = models.load_ssp_cube("/home/vital/Downloads/pySB99_SESAMME_Cube_v4/pySB99_SSP_Grid_v4.fits")
# iontable = models.load_ionization_table("/Users/sveash/Library/CloudStorage/Dropbox/SESAMME_HIs/SESAMME/pySB99_SESAMME_SSPs_v1/pySB99_ionization_Qtable_v1.txt")

# Natively comes in flux density units (erg s-1 cm-2 A-1)...
specfile = Table.read(spec_pname, delimiter=' ', format='ascii')
# wl = specfile['WL']
# flux = specfile['FLUX']
# flux_err = specfile['ERROR']





# But rescaling to luminosity density units (L_Sun A-1) allows SESAMME to infer a stellar mass for the cluster
# lum = flux * 4*np.pi * (distance*u.Mpc.to(u.cm))**2 / 3.83e33
# lum_err = flux_err * 4*np.pi * (distance*u.Mpc.to(u.cm))**2 / 3.83e33

# units of Vital's grids are erg s^-1 AA^-1
# lum = flux * 4 * np.pi * (distance * u.Mpc.to(u.cm)) ** 2
# lum_err = flux_err * 4 * np.pi * (distance * u.Mpc.to(u.cm)) ** 2

# ---set the masks to be applied to data
# windowlist = np.array([[np.min(wl),1000],[1010,1147],[1163,1225],[1250,1305],[1325,1528],[1600,1690],[1748, np.max(wl)] ])
# windowlist = np.array(
#     [[np.min(wl), 1101.61], [1104.1, 1119.01], [1120.5, 1135.1], [1137.99, 1146], [1150, 1158.7], [1162.35, 1178],
#      [1184.7, 1222.97], [1245.9, 1266], [1289, 1308], [1317, 1344.76], [1354, 1367], [1381.2, 1404],
#      [1422.21, np.max(wl)]])
# # 1114],[1116.9
# mask = models.get_mask(windowlist, wl)

# Set the initial positions of the walker ensemble
# logAge= 6.74
# logZ = -2.398
# ebv= 0.04
# logA = -0.73
initial = [6.74, -2.398, 0.04, -0.73]

# -- this is the likelihood function that is being minimized
nll = lambda *args: -log_likelihood_model(*args)

solution = minimize(nll, initial, method='Nelder-Mead', bounds=((6, 7.),
                                                                (np.log10(0.001), np.log10(0.004)),
                                                                (0.0, 0.2),
                                                                (-1.5, 2)),
                    args=(wl, lum, lum_err, mask, modelcube, None, False),
                    options={'maxiter': 8000})

initial_optimized = solution.x.tolist()

print('initial_optimized', initial_optimized)
# print(initial_optimized)

# Set the dimensionality of the emcee walker ensemble to N x 4
stats.set_walker_size(128)

# Set the desired chain length
stats.set_chain_size(2000)

# Set the initial positions of the walker ensemble
# here we are using the values inferred from the minimization above
stats.set_initial_positions(initial_optimized)
# print(stats.initial_pos)

# prior_lowbounds = [6.0, np.log10(0.001), 0.0, -1.5]
# prior_highbounds = [7., np.log10(0.01), 0.6, 2]
prior_lowbounds = [6.0, np.log10(0.001), 0.0, -1.5]
prior_highbounds = [7., np.log10(0.01), 0.6, 2]

stats.set_prior_bounds(stats.prior_dict, prior_lowbounds, prior_highbounds)


obj = 'NGC4861'
filename = "NGC4861_pySB99.h5"
runname = 'Run_2000_v2'

#--- Note: set the variable below to False if using SB99 models, and True if you are using BPASS.
#--- This variable dictates if SESAMME should include or exclude nebular continuum in the models.
neb_continuum = False
stats.run_sesamme(filename, runname, wl, lum, lum_err, modelcube, None, mask, neb_continuum)
# [7.0, -2.3979400086720375, 0.2, -1.5]
#
print(prior_lowbounds)
#
# reader = emcee.backends.HDFBackend(filename, name = runname)
# samples = reader.get_chain()

# Save Statistics
print(f'Loading {filename}:')
reader = emcee.backends.HDFBackend(filename, name=runname)
tau = np.array([121.55176189, 204.45805254,  92.25215704, 148.32895022])
thin_num = int(tau.max()/2.)
flat_samples = reader.get_chain(discard=40, thin=thin_num, flat=True)
grid_samples = samples = reader.get_chain(discard=40, thin=thin_num, )

vis.print_stats(flat_samples)

fig, axes = plt.subplots(4, figsize=(9, 7), sharex=True)

labels = ["log(age/yr)", r"log(Z/Z$_{\odot}$)", "E(B-V)", "log(A)"]

for i in range(stats.ndim):
    ax = axes[i]
    ax.plot(grid_samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(grid_samples))

    ax.set_ylabel(labels[i])
    ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number")
plt.show()

vis.save_stats(flat_samples, output_path='./', run_name=f'ngc_{runname}')
fig = corner.corner(flat_samples, **CORNER_CFG)
fig.savefig(f'{obj}_{runname}_corner.png', dpi=300)

vis.plot_samples(wl, lum, get_masked_bands(mask, wl), flat_samples, add_nebular=False)
