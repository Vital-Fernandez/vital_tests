import re
import lime
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from pyneb import RedCorr
from astropy.convolution import Gaussian1DKernel
from specutils import Spectrum1D, Spectrum
from specutils.manipulation.smoothing import convolution_smooth
import astropy.units as u
from astro.stsci.tools import deredden
from astropy.nddata import StdDevUncertainty
from pyneb import RedCorr


def extinction_corr(R_V, A_V, in_wave, in_flux, in_err, obj, plot_results=False):
    E_BV = A_V / R_V

    # UV correction flux
    flux_arr_gordon, err_arr_gordon = deredden(R_V, E_BV, in_wave, in_flux, in_err)

    if plot_results:
        rc = RedCorr(R_V=R_V, law='G03 LMC', E_BV=E_BV)
        red_corr = rc.getCorr(spec.wave)

        fig, ax = plt.subplots()
        ax.step(in_wave, in_flux, where='mid', label='Observed')
        ax.step(in_wave, in_flux * red_corr, where='mid', label='Gordon 2003 (Optical)', linestyle='--')
        ax.step(in_wave, flux_arr_gordon, where='mid', label='Gordon 2009 (UV)', linestyle=':')
        ax.set_title(f'{obj} MW reddening correction E(B-V) = {E_BV:0.3f}', fontsize=18)
        ax.legend(fontsize=16)
        plt.show()

    return flux_arr_gordon, err_arr_gordon


def smooth_corr(smooth_sigma, in_wave, in_flux, in_err, obj, plot_results=False):

    kernel = Gaussian1DKernel(stddev=smooth_sigma)

    spec1d = Spectrum(spectral_axis=wave_arr * u.AA, flux=in_flux * u.erg / (u.cm**2 * u.s * u.AA),
                      uncertainty=StdDevUncertainty(in_err * u.erg / (u.cm**2 * u.s * u.AA)))
    spec1d_s = convolution_smooth(spec1d, kernel)

    StdDevUncertainty()
    if plot_results:
        label = 'Smoothed spectrum'
        fig, ax = plt.subplots()
        ax.step(in_wave, in_flux, where='mid', label='Input spectrum')
        ax.step(in_wave, spec1d_s.flux, where='mid', label=label, linestyle='--')
        ax.set_title(f'{obj} smoothing MODEL(SB99) - COS  ($\sigma$={smooth_sigma:0.3f})', fontsize=18)
        ax.legend(fontsize=16)
        plt.show()

    return spec1d_s.flux.value, spec1d_s.uncertainty.array


def rebin_corr(input_wave, input_flux, input_err, ssp_wave, obj, plot_results=False):
    """
    Returns the slice of arr2 that spans the same range as arr1,
    using searchsorted for efficiency on sorted arrays.
    """

    # Get the interval for interpolation
    idx_start, idx_end = np.searchsorted(ssp_wave, (input_wave[0], input_wave[-1]))
    disp_intvl = ssp_wave[idx_start:idx_end]

    # New
    spec = lime.Spectrum(wave_arr, input_flux, input_err, redshift=0)
    wave_b, flux_b, err_b = spec.retrieve.rebinned(disp_intvl=disp_intvl)

    if plot_results:
        fig, ax = plt.subplots()
        ax.step(input_wave, input_flux, where='mid', label='Input spectrum')
        ax.step(wave_b, flux_b, where='mid', label='Rebineed spectrum', linestyle='--')
        ax.set_title(f'{obj} rebining to pySB99 resolution', fontsize=18)
        ax.legend()
        plt.show()

    return wave_b, flux_b, err_b


def lumin_corr(input_flux, input_err, distance_mpc):

    distance_cm = (distance_mpc * u.Mpc).to(u.cm).value
    conversion = 4 * np.pi * distance_cm**2 / 3.83e33
    lum = input_flux * conversion
    lum_err = input_err * conversion

    return lum, lum_err


# Data location
obs_folder = Path('/home/vital/Astrodata/STScI/LyC_leakers_COS')
project_folder = Path('/home/vital/Dropbox/Astrophysics/Data/STScI_projects')
ssp_wave_arr = np.loadtxt('/home/vital/Dropbox/Astrophysics/Data/STScI_projects/LyC_leakers_COS/SSPs/pyStarburst99_wave_arr.txt')
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
galaxy_distances = cfg_sample['galaxy_distances']


# Loop through the targets and generate the HASP files
for i, obj in enumerate(target_list):

    if i > 0: #11 13

        # Object folders the files
        spec_folder = project_folder/'LyC_leakers_COS'/'metals_line_fitting'
        stellar_folder = project_folder/'LyC_leakers_COS'/'SSPs'
        spec_pname = spec_folder/f'{obj}_rebinned_LymanAlphaCorrected.txt'
        mask_pname = stellar_folder/'stellar_masks'/f'{obj}_mask_intervals.txt'
        print(f'{i}) Galaxy {obj}')

        # Load the spectrum
        spec = lime.Spectrum.from_file(spec_pname, instrument='text')
        wave_arr = spec.wave_rest.data
        flux_arr = spec.flux.data * spec.norm_flux * (1 + spec.redshift)
        err_arr = spec.err_flux.data * spec.norm_flux * (1 + spec.redshift)

        # Get the extinction
        flux_dered_arr, err_dered_arr = extinction_corr(R_V=3.4, A_V=cfg_sample['galaxy_extinction'][obj],
                                                        in_wave=wave_arr, in_flux=flux_arr, in_err=err_arr, obj=obj,
                                                        plot_results=False)

        # Smooth the spectrum
        smooth_sigma = np.sqrt(np.square(0.4 / 2.355) - np.square(0.06 / 2.355))
        # smooth_sigma = np.sqrt(np.square(cfg_sample['LyC_acq_image_fwhm_pixels'][obj] / 2.355) - np.square(0.06 / 2.355))
        flux_s_arr, err_s_arr = smooth_corr(smooth_sigma, in_wave=wave_arr, in_flux=flux_dered_arr, in_err=err_dered_arr, obj=obj, plot_results=False)
        flux_s_arr, err_s_arr = flux_dered_arr, err_dered_arr

        # Mask the spectrum
        spec = lime.Spectrum(wave_arr, flux_s_arr, err_s_arr, redshift=0)
        wave_rebin_arr, flux_rebin_arr, err_rebin_arr = rebin_corr(wave_arr, flux_s_arr, err_s_arr, ssp_wave_arr, obj=obj, plot_results=False)

        # Convert to luminosity
        flux_lum, err_lum = lumin_corr(flux_rebin_arr, err_rebin_arr, distance_mpc=galaxy_distances[obj])

        # Mask the spectrum
        spec = lime.Spectrum(wave_rebin_arr, flux_lum, err_lum, redshift=0, units_flux='L_sun')
        obj_mask_spec = spec.retrieve.spectrum(mask_intvls=mask_pname, norm_flux=1)
        # obj_mask_spec.plot.spectrum()

        # Save the spectrum
        obj_mask_spec.save_spectrum(stellar_folder/f'{obj}_pyStarburst99_spectrum_Lsun.txt')

