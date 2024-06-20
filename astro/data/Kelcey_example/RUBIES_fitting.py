import numpy as np
import lime
from astropy.io import fits
from pathlib import Path

file_address = 'jw04233005002_CEERS-FULL-V2_CLEAR_PRISM_s08488_x1d.fits'
redshift = 8.69

# Spectrum from .fits file (MJy and micrometers)
spec = lime.Spectrum.from_file(file_address, 'nirspec', redshift=redshift)
spec.unit_conversion('AA', 'FLAM', norm_flux=1e-15)
spec.plot.spectrum(rest_frame=True)

# Opening .fits file
with fits.open(file_address) as hdu1:
    datap1 = hdu1[1].data
    wave_1d = datap1['WAVELENGTH']
    flx_1d = datap1['FLUX']
    flx_es_1d = datap1['FLUX_ERROR']

# Adding nan entries for lower limit
emax = 1e-12
idcs_nan = flx_1d < emax
flx_1d[idcs_nan] = np.nan
flx_es_1d[idcs_nan] = np.nan
pixel_mask = np.isnan(flx_1d)

spec_v2 = lime.Spectrum(wave_1d, flx_1d, flx_es_1d, redshift=redshift, units_wave='um', units_flux='MJy', pixel_mask=pixel_mask)
spec_v2.unit_conversion('AA', 'FLAM', norm_flux=1e-15)
spec_v2.plot.spectrum(rest_frame=True)

# Comparing maximum flux in both files
for i, obj in enumerate([spec, spec_v2]):
    idx_max = np.nanargmax(obj.flux.data)
    print(f'Opening ({i}) Coordiante max: {idx_max} , max flux {obj.flux[idx_max]}, max err {obj.err_flux[idx_max]}')

# Prepare a bands file
bands_file = '8488_bands.txt'
if not Path(bands_file).is_file():
    bands_clean_df = lime.line_bands(wave_intvl=spec, vacuum=True)
    lime.save_frame(bands_file, bands_clean_df)

# Adjust the bands:
spec.check.bands(bands_file, ref_bands=lime.line_bands(wave_intvl=spec, vacuum=True))

# Fit the lines
fit_cfg_file = 'RUBIES_cfg.toml'
spec.fit.frame(bands_file, fit_cfg_file, id_conf_prefix='8488')
spec.plot.spectrum(rest_frame=True)

# Save the measurements
spec.save_frame('8488_fluxes_table.txt')
spec.save_frame('rubies_fluxes.fits', page='8488')
