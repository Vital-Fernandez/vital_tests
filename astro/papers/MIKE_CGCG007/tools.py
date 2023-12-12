import numpy as np
from astropy.io import fits
import lime

def read_mike_fits(file_path):

    with fits.open(file_path) as hdul:
        wave_m, flux_m, err_m = hdul[1].data, hdul[2].data, hdul[3].data

    return wave_m, flux_m, err_m


def mike_load_function(log_df, id_spec, **kwargs):

    file_path = log_df.loc[id_spec, 'file_path']
    order = log_df.loc[id_spec, 'order']
    z_obj = kwargs['redshift']

    wave_m, flux_m, err_m = read_mike_fits(file_path)

    wave_array, flux_array, err_array = wave_m[order, :], flux_m[order, :], err_m[order, :]
    spec_j = lime.Spectrum(wave_array, flux_array, err_array, redshift=z_obj)

    return spec_j


def nirspec_load_function(log_df, id_spec, **kwargs):

    z_obj = log_df.loc[id_spec].redshift
    file_spec = Path(kwargs['fits_folder'])/log_df.loc[id_spec].file_path

    # 1d files
    if file_spec.as_posix().endswith('_x1d_masked.fits'):
        wave, flux, err, header = load_nirspec_fits(file_spec)
        norm_flux = kwargs['norm_flux']

        mask = np.isnan(err) & np.isnan(flux)
        objSpec = lime.Spectrum(wave, flux, err, redshift=z_obj, units_wave='um', units_flux='Jy', pixel_mask=mask)
        objSpec.unit_conversion(units_wave='A', units_flux='Flam', norm_flux=norm_flux)
        objSpec.header = header

    # 2d files
    else:
        wave, flux, err, header = load_nirspec_fits(file_spec)
        objSpec = wave, flux, err, header

    return objSpec