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
# from specsy.models.stellar_synthesis import get_initial_values, build_grid, set_prior_bounds, build_model
from specsy.models.stellar_synthesis import SSP_sampler, load_fit_results, plot_best_fit, load_ionization_frame, summary_ssp
import arviz as az

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
ssp_wave_arr = np.loadtxt('/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs/pyStarburst99_wave_arr.txt')
stellar_folder = project_folder / 'LyC_leakers_COS' / 'SSPs'
SpecSy_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS/SpecSy_fittings')

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
ssp_cube.plot_age_metallicity()

# Boundary conditions
pyStaburst99_bounds = np.array([(6, 7.5), np.log10([0.00001, 0.02]), (0.0, 0.5),  (-2, 20)])

# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    if i == 6: #11 13

        # Declare inputs
        print(f'{i}) Galaxy {obj}')
        label_fit = '_1population_pyStarburst_v1'
        input_pname = stellar_folder / f'{obj}_pyStarburst99_spectrum_Lsun.txt'
        output_trace = SpecSy_folder / f'{obj}_{label_fit}_trace'

        # Load the spectrum
        spec = lime.Spectrum.from_file(input_pname, instrument='text', norm_flux=1)
        spec.plot.spectrum(ax_cfg = {'title': f'Galaxy: {obj}'}, rest_frame=True)
        wave_arr, flux_arr, err_arr, mask_arr = spec.retrieve.spectrum(return_arrays=True)
        mask_arr = ~mask_arr

        # Get the object SSP
        obj_ssp = ssp_cube.to_stellar_binaries(disp_intvl=wave_arr)

        # Grid nodes on the observed wavelengths (one get_spectrum call per node), plus the nebular continuum
        sampler = SSP_sampler(spec, obj_ssp, red_law='Calzetti', r_v=3.1)
        sampler.prepare_inputs(wave_arr, flux_arr, err_arr, mask_arr, prior_bounds=pyStaburst99_bounds, add_nebular=False, ion_frame=None)
        trace = sampler.sample(draws=100, tune=1000, spread_step=True, progressbar=True)
        sampler.summary()

        # Save the results
        sampler.save_trace(output_trace)

        # Load the results
        trace = az.from_netcdf(output_trace)

        # Load the results
        trace_data = az.from_netcdf(output_trace)
        results_df = summary_ssp(trace_data)
        print(results_df)

        # Plot the results
        plot_best_fit(trace_data)