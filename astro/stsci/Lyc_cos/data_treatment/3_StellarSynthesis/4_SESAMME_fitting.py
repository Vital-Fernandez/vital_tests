import numpy as np
import re
import emcee
import lime
from pathlib import Path
from specsy.models.ssp import StellarBinaries
import sesamme.mcmc as stats
from sesamme import models, vis
from sesamme.models import load_ionization_table, load_ssp_cube
from scipy.optimize import minimize
import corner
from matplotlib import pyplot as plt
from astropy.io import fits
from astro.stsci.tools import get_masked_bands, CORNER_CFG, plot_chains, plot_sesamme_samples


def log_likelihood_model(theta, wl, y, yerr, mask, model_cube, ion_table, add_nebular):
    y_model = models.get_model(theta, wl, model_cube, ion_table, add_nebular)

    masked_spec = np.array((y[mask] - y_model[mask]) / yerr[mask])
    masked_err = np.array(np.sqrt(2) * np.sqrt(np.pi) * yerr[mask])

    resid = -0.5 * (np.dot(masked_spec, masked_spec) + np.log(np.dot(masked_err, masked_err)))

    return resid


def get_initial_values(initial, wl, lum, lum_err, mask, modelcube, bounds, show_solution=True):

    # -- this is the likelihood function that is being minimized
    nll = lambda *args: -log_likelihood_model(*args)

    solution = minimize(nll, initial, method='Nelder-Mead', bounds=bounds,
                        args=(wl, lum, lum_err, mask, modelcube, None, False),
                        options={'maxiter': 8000})

    if show_solution:
        print(solution)

    initial_optimized = solution.x.tolist()

    return initial_optimized


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
ssp_wave_arr = np.loadtxt('/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs/pyStarburst99_wave_arr.txt')
stellar_folder = project_folder / 'LyC_leakers_COS' / 'SSPs'
SESAMME_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/SESAMME_fittings')

# Cfg file
cfg_sample = lime.load_cfg('../../Lyc_cos.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))
sample_df = sample_df.loc[~sample_df["filepath"].str.contains(pattern, na=False)]

# Order for objects
target_list = sample_df.object.unique()
target_list = ['SBS0335052', 'IZw18', 'SBS1415437', 'SBS1159545', 'UM461', 'Pox186',
               'UGCA281', 'NGC1705',  'NGC4861',  'Haro2', 'He2-10', 'MRK1450',         # Codo
               'UGC4483', 'VIIZw403',  'NGC2366',                                       # Dent
               'Haro11_A', 'Haro11_B', 'Haro11_C']

# Load the cube
# ssp_cube_pname = Path('/home/vital/Downloads/pySB99 _SESAMME_Cube_v4/pySB99_SSP_Grid_v4.fits')
ssp_cube_pname = Path('/home/vital/Downloads/pySB99_SESAMME_Cube_v5/pySB99_SSP_Grid_v5.fits')
ssp_cube = StellarBinaries.from_fits(ssp_cube_pname, source='pystarburst99')
# ssp_cube.plot_age_metallicity()

# Boundary conditions
# pyStaburst99_bounds = ((6, 7.), (np.log10(0.001), np.log10(0.004)), (0.0, 0.2),  (-1.5, 2))
# pyStaburst99_bounds = ((6, 7.), np.log10([0.0004, 0.02]), (0.0, 0.2),  (-1.5, 2))
pyStaburst99_bounds = np.array([(6, 7.5), np.log10([0.00001, 0.02]), (0.0, 0.5),  (-1.5, 2)])

# Extinction law
models.set_ext_law('Calzetti')


# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    if i == 6: #11 13

        # Declare inputs
        print(f'{i}) Galaxy {obj}')
        input_pname = stellar_folder / f'{obj}_pyStarburst99_spectrum_Lsun.txt'
        spec = lime.Spectrum.from_file(input_pname, instrument='text')
        spec.plot.spectrum(ax_cfg = {'title': f'Galaxy: {obj}'}, rest_frame=True)

        wave_arr = spec.wave_rest.data
        flux_arr = spec.flux.data * spec.norm_flux
        err_arr = spec.err_flux.data * spec.norm_flux
        mask = ~spec.wave_rest.mask

        # Prepare stellar inputs
        runname = f'v4_default_pystarburst99'
        fname_sesamme = SESAMME_folder/runname/f'{obj}_SESAMME.h5'
        fname_sesamme.parent.mkdir(parents=True, exist_ok=True)
        obj_ssp_pname = SESAMME_folder/runname/f'{obj}_object_ssp_{ssp_cube.source}.fits'
        if not obj_ssp_pname.is_file():
            ssp_cube.to_fits(obj_ssp_pname, disp_intvl=wave_arr)

        obj_cube = models.load_ssp_cube(obj_ssp_pname)

        # Fit the initial values
        p0_opt = cfg_sample['p0_SESAMME'].get(obj)
        if p0_opt is None:
            p0_opt = [6.5, np.log10(0.001), 0.05, -0.7]
            p0_opt = get_initial_values(p0_opt, wave_arr, flux_arr, err_arr, mask, obj_cube, pyStaburst99_bounds,
                                        show_solution=True)
            print(f'Computed Optimized inputs: {p0_opt}')
        else:
            print(f'User input values: {p0_opt}')

        # Set the sampler configuration
        stats.set_walker_size(128)
        stats.set_chain_size(4000)
        stats.set_initial_positions(p0_opt)

        prior_lowbounds = pyStaburst99_bounds[:, 0]
        prior_highbounds = pyStaburst99_bounds[:, 1]
        stats.set_prior_bounds(stats.prior_dict, prior_lowbounds, prior_highbounds)

        # Run the sampler
        stats.run_batch_sesamme(fname_sesamme, runname, wave_arr, flux_arr, err_arr, obj_cube, ion_table=None, mask=mask, add_nebular=False)

        # Load the results
        reader = emcee.backends.HDFBackend(fname_sesamme, name=runname)
        flat_samples = reader.get_chain(discard=2000, thin=10, flat=True)
        grid_samples = samples = reader.get_chain(discard=2000, thin=10)

        # Plot the results
        fig = corner.corner(flat_samples, **CORNER_CFG)
        fig.savefig(SESAMME_folder/runname/f'{obj}_{runname}_corner.png', dpi=300)

        plot_sesamme_samples(wave_arr, flux_arr, get_masked_bands(~mask, wave_arr), flat_samples, add_nebular=False,
                             model_cube = obj_cube, savefile_name=SESAMME_folder/runname/f'{obj}_{runname}_fitted_spec.pdf')

        vis.save_stats(flat_samples, output_path=SESAMME_folder / runname, run_name=f'{obj}_{runname}')
        plot_chains(grid_samples, stats, fname=SESAMME_folder/runname/f'{obj}_{runname}_chains_spec.pdf')

