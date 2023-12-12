import numpy as np
from pathlib import Path
from mpdaf.obj import Cube
from astropy.wcs import WCS
import lime


def read_muse_cube(file_address):

    cube_obj = Cube(filename=str(file_address))
    header = cube_obj.data_header

    dw = header['CDELT3']
    w_min = header['CRVAL3']
    nPixels = header['NAXIS3']
    w_max = w_min + dw * nPixels
    wave_array = np.linspace(w_min, w_max, nPixels, endpoint=False)

    return wave_array, cube_obj, header


# Inputs
cfg_file = 'NGC1386_muse.toml'
cube_file = Path('D:/Dropbox/Burocracy/Trabajo/2023_Schmidt_Postdoc/Proposal/subcube_1386.fits')
bands_file = 'NGC1386_bands.txt'
log_file = 'NGC1386_log.txt'

# Load configuration
cfg = lime.load_cfg(cfg_file)
norm_flux = 1
z_obj = 0.00259

# Load cube
wave_array, cube, hdr = read_muse_cube(cube_file)
flux_cube = cube.data.data * norm_flux
mask_pixel_cube = np.isnan(flux_cube)
wcs = WCS(hdr)

# Create MUSE
NGC1386 = lime.Cube(wave_array, flux_cube, redshift=z_obj, norm_flux=norm_flux,
                    wcs=wcs, pixel_mask=mask_pixel_cube)

# NGC1386.check.cube('H1_6563A', rest_frame=True)

spec_BH = NGC1386.get_spectrum(23, 25)
# spec_BH.plot.spectrum(rest_frame=False)

# Create the bands
# spec_BH.check.bands('NGC1386_bands.txt')

target_lines = ['O3_5007A_b', "H1_6563A_b"]
spec_BH.fit.frame(bands_file, cfg, line_list=target_lines, plot_fit=False, default_conf_prefix='no_decimal')
spec_BH.save_log(log_file)


fig_cfg = {'figure.dpi': 400, 'figure.figsize': (10, 8)}
spec_BH.plot.bands('H1_6563A_b', output_address='Halpha_fitting.png', fig_cfg=fig_cfg)