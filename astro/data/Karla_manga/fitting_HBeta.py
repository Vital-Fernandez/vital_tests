import numpy as np
import lime
from pathlib import Path
from astropy.io import fits

path_fits = Path('NGC5471_C_34-28-LR-B.fits')
configFile = Path('ejemplo_ajuste.cfg')
bands_file = Path('MEGARA_bands.txt')
obs_cfg = lime.load_cfg(configFile)
z_obj = 0.00091

# Load the spectrum
frame_idx = 0
with fits.open(path_fits) as hdul:
    flux, hdr = hdul[frame_idx].data, hdul[frame_idx].header
    w_min = hdr['CRVAL1']
    dw = hdr['CDELT1']  # dw (Wavelength interval per pixel)
    pixels = hdr['NAXIS1']  # nw number of output pixels
    w_max = w_min + dw * pixels
    wave = np.linspace(w_min, w_max, pixels, endpoint=False)

gp_spec = lime.Spectrum(wave, flux, redshift=z_obj, norm_flux=1, units_flux='Jy')
gp_spec.convert_units(units_flux='Flam', norm_flux=1e-17)
gp_spec.plot.spectrum(rest_frame=True)

gp_spec.fit.frame(bands_file, fit_conf=obs_cfg['LR-B_line_fitting'], lines_list=['O3_4959A_b'], line_detection=False)
gp_spec.save_log(f'LR-B_log.txt')

fig_conf = {'figure.dpi': 200,  'lines.linewidth': 0.5, 'font.size': 8}
gp_spec.plot.band(rest_frame=True, fig_cfg=fig_conf)


# gp_spec.save_log(f'LR-B_log.txt')


# line = 'N2_5755A_b'
# band_edges = np.array([5700.000000, 5720.000000, 5730.000000, 5785.000000, 5800.000000, 5820.000000])
# gp_spec.fit.band(line, band_edges, fit_conf=obs_cfg['LR-V_line_fitting'])
# gp_spec.plot.band()
# gp_spec.save_log(f'N2_5755A_log.txt')

# line = 'H1_6562.7192A_b'
# band_edges = np.array([6438.03, 6508.66, 6525.10, 6616.0, 6627.70, 6661.82])
# gp_spec.fit.band(line, band_edges, fit_conf=obs_cfg['LR-V_line_fitting'])
# gp_spec.plot.band()

# gp_spec.save_log(f'manga_2compsHR_plus_wide_line_fitting.txt')

# line = 'H1_6562.7192A_b'
# gp_spec.fit.band(line, band_edges, fit_conf=obs_cfg['manga_3compsHR_line_fitting'])
# gp_spec.plot.band()
#
# line = 'H1_6562.7192A_b'
# gp_spec.fit.band(line, band_edges, fit_conf=obs_cfg['manga_3compsHR_plusOne_line_fitting'])
# gp_spec.plot.band()

#
# print(lime.__version__)