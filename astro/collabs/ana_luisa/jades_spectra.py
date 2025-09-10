import numpy as np
from astropy.io import fits
import lime
from tests.test_tools import redshift


def read_rubies(fname, z_obj=None):
    with fits.open(fname) as spectrum_data:
        data = spectrum_data[1].data
        wave = np.array(data['WAVELENGTH'])
        flux = np.array(data['FLUX'])
        flux_err = np.array(data['FLUX_ERR'])

        spec_lime = lime.Spectrum(wave, flux, flux_err, units_wave='um', units_flux='FLAM', redshift=z_obj)
        spec_lime.unit_conversion(wave_units_out='AA', flux_units_out='FLAM')

    return spec_lime

fname1 = 'hlsp_jades_jwst_nirspec_goods-s-deephst-00004282_clear-prism_v1.0_x1d.fits'
fname2 = 'hlsp_jades_jwst_nirspec_goods-s-deephst-00004282_f290lp-g395h_v1.0_x1d.fits'

# Open the spectra
spec1 = read_rubies(fname1, z_obj=3.611)

# # Generate the bands
# bands = spec1.retrieve.line_bands(vacuum_waves=True, band_vsigma=250)Z
# lime.save_frame('example_bands.txt', bands)

spec1.fit.bands('H1_6565A', bands='example_bands.txt')
spec1.plot.bands(rest_frame=True, output_address='Halpha.png')
spec1.plot.spectrum()
