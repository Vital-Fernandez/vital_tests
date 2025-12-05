import re
from pathlib import Path
import numpy as np
import lime
from astropy.io import fits
from matplotlib import pyplot as plt, rc_context

# Data location
obs_folder = Path('/home/vital/Astrodata/STScI')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')

# Cfg file
cfg_sample = lime.load_cfg(project_folder/'samples.toml')

# Sample file
sample_df = lime.load_frame(project_folder/'stsci_samples_v1.csv', levels=['sample', 'id', 'offset_id', 'state'])
pattern = "|".join(map(re.escape, sum(cfg_sample['excluded_files'].values(), [])))

idcs_manual = sample_df.index.get_level_values('state') == 'aspec_manual'
obj_list = sample_df.loc[idcs_manual, 'object'].unique()

for i, obj_name in enumerate(obj_list):

    # Get subspec
    idcs_cspec = (sample_df.object == obj_name) & (sample_df.index.get_level_values('state') == 'cspec_manual')
    sub_sample = sample_df.loc[idcs_cspec]

    # Manual cspec
    if sub_sample.index.size == 0:
        print(f'-{obj_name} NO SPECTRA')

    else:
        redshift = cfg_sample['Galaxy_redshifts'][obj_name]
        file_path = obs_folder / sample_df.loc[idcs_manual & (sample_df.object == obj_name), 'filepath'].values[0]
        print(f'- {obj_name} at z = {redshift} ({file_path.stem})')

        with rc_context({'figure.figsize':(6, 4), 'figure.dpi': 150}):

            cmap = plt.get_cmap("plasma_r")
            color_list = cmap(np.linspace(0, 1, len(sub_sample)))

            spec = lime.Spectrum.from_file(file_path, instrument='COS', redshift=redshift)
            wave, flux, err = spec.retrieve.rebinned(pixel_number=6, constant_pixel_width=True)

            fig, ax = plt.subplots()
            ax.step(spec.wave, spec.flux, where='mid', color='black', label=f'{file_path.stem}')
            ax.step(wave, flux, where='mid', color='red', label=f'Rebinned', linestyle='--')

            # for j, idx in enumerate(sub_sample.index):
            #     file_path = obs_folder / sub_sample.loc[idx, 'filepath']
            #     c_spec = lime.Spectrum.from_file(file_path, instrument='COS', redshift=redshift)
            #     wave_i, flux_i, err_i = c_spec.retrieve.rebinned(pixel_number=6, constant_pixel_width=True)
            #     ax.step(wave_i, flux_i, where='mid', label=f'{file_path.stem}', color=color_list[j], linestyle='--')

            ax.update({'title': f'{obj_name}'})
            ax.legend()
            plt.tight_layout()
            plt.show()


import lime
# import numpy as np
# from scipy import stats
# import spectres
#
#
# def bin_data(wave, flux, err, opt_elem):
#
#     '''
#
#     This script bins data
#     by 1 resolution element for COS (6 pixes).
#
#     # get errors on binned flux #
#     1) loop through the bin numbers,
#     2) find the index in the binnumber array,
#     3) square the error (for the corresponding index)
#     4) sum the squared errors
#     5) get the square root of the summed errors
#     6) append to binn_err array
#
#     '''
#
#     # This is in Angstroms per resel
#     if opt_elem == 'G130M':
#         bin_width =  0.05982
#     elif opt_elem == 'G160M':
#         bin_width = 0.07338
#
#     # Bin data
#     bin_edges = np.arange(wave[0],wave[-1], bin_width)
#     binned_flux, edges, binnumber = stats.binned_statistic(wave, flux, statistic='mean', bins=bin_edges)
#     binned_wave = bin_edges[:-1]+(bin_width/2.)
#
#     # get unique bin numbers #
#     uni_bin = np.unique(binnumber)
#     err_binned = []
#
#     for binnum in uni_bin[:-1]:
#         index_bin = np.where(binnumber==binnum)
#         errors_bin = err[index_bin]
#         err_bin = np.sqrt(np.sum(errors_bin**2))/(len(errors_bin))
#         err_binned.append(err_bin)
#
#     return binned_wave, binned_flux, err_binned
#
# file_path = '/home/vital/Astrodata/STScI/LyC_leakers_COS/Direct_downloads/LF9G01010/hst_17515_cos_mrk-209_g130m_lf9g01_cspec.fits'
# spec = lime.Spectrum.from_file(file_path, instrument='cos', redshift=0.000932, norm_flux=1e-17)
#
# spec.plot.spectrum(show_err=True, in_fig=None)
#
# wave_svea, flux_svea, err_svea = bin_data(spec.wave, spec.flux, spec.err_flux, opt_elem='G130M')
# flux_carnall, err_carnall = spectres.spectres(wave_svea, spec.wave, spec.flux, spec.err_flux)
#
# spec.plot.ax.step(wave_svea, flux_svea, label='Svea rebin', where='mid')
# spec.plot.ax.fill_between(x=wave_svea, y1=(flux_svea - err_svea), y2=(flux_svea + err_svea), step='mid', alpha=0.5, color='blue', ec=None)
#
# spec.plot.ax.step(wave_svea, flux_carnall, label='Carnall rebin', where='mid')
# spec.plot.ax.fill_between(x=wave_svea, y1=(flux_carnall - err_carnall), y2=(flux_carnall + err_carnall), step='mid', alpha=0.5, color='orange', ec=None)
#
# spec.plot.ax.legend()
# spec.plot.show()