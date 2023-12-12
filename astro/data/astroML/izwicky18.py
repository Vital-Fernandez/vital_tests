import lime
from pathlib import Path
from astropy.io import fits
import numpy as np

def import_fits_data(file_address, frame_idx):

    # Open fits file
    with fits.open(file_address) as hdul:
        data, header = hdul[frame_idx].data, hdul[frame_idx].header

    # Check instrument
    if 'INSTRUME' in header:
        if 'ISIS' in header['INSTRUME']:
            instrument = 'ISIS'

    # William Herschel Telescope ISIS instrument
    if instrument == 'ISIS':
        w_min = header['CRVAL1']
        dw = header['CD1_1']  # dw = 0.862936 INDEF (Wavelength interval per pixel)
        pixels = header['NAXIS1']  # nw = 3801 number of output pixels
        w_max = w_min + dw * pixels
        wave = np.linspace(w_min, w_max, pixels, endpoint=False)

    return wave, data, header


file_address = Path(f'D:\Dropbox\Astrophysics\Tools\LineMesurer\IZwicky18\IZW18_Blue_cr_f_t_w_e__narrowSlits_fglobal.fits')

wave, flux, hdr = import_fits_data(file_address, 0)

spec = lime.Spectrum(wave, flux[1, :], redshift=0.00095)
output_folder = Path(f'D:/Dropbox/Astrophysics/Seminars/Umich_introduction_2023')

spec.plot.spectrum()   #spectrum(output_address=output_folder/'Izwicky18.png')
# spec.plot.bands(6312, output_address=output_folder/f'line_redshift_{6312}.png')
# spec.fit.bands(6312)
# spec.plot.bands(output_address=output_folder/f'line_redshift_{6312}_fit.png')