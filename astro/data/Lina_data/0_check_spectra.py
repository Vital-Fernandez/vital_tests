from astropy.io import fits
from pathlib import Path
import lime

def read_mike_fits(file_path):

    with fits.open(file_path) as hdul:
        wave_m, flux_m, err_m = hdul[1].data, hdul[2].data, hdul[3].data

    return wave_m, flux_m, err_m


file_addres = f'J082652+182052.fits'
bands_address = f'J082652+182052_bands.txt'
config_address = f'fit_conf.toml'

fits.info(file_addres)
data = fits.getdata(file_addres, ext=0)
hdr = fits.getheader(file_addres, ext=0)
wave, flux = data[0], data[1]

spec = lime.Spectrum(wave, flux)
spec.fit.bands('O3-3comps_5007A_b', bands_address, config_address)
# spec.fit.bands('O3-exp_5007A_b', bands_address, config_address)
spec.plot.bands()
