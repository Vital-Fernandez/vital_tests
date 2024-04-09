import numpy as np
import lime
from pathlib import Path

spec_address = './Hf2_2.fits'
mask_plot = './Hf2_2_mask.png'
cfgFile = './Hf2_2.toml'
spatial_mask = './Hf2_2_mask.fits'

# Cargamos las lineas para leer
spatial_mask_file = Path(spatial_mask)
bands_file = Path('./Lineas_muse_minimo.txt')
bands = lime.load_frame(bands_file)
obs_cfg = lime.load_cfg(cfgFile)

output_lines_log_file = Path('./Hf2_2_log.fits')

# Observation properties
z_obj = obs_cfg['sample_data']['redshift']
norm_flux = obs_cfg['sample_data']['norm_flux']

hf2 = lime.Cube.from_file(spec_address, instrument='MUSE', redshift=z_obj)

# spec_malo = hf2.get_spectrum(189, 322)
# spec_malo.plot.spectrum()
#
# matched_bands = spec_malo.line_detection(bands=bands_file)
# print(matched_bands)
#
# spec_malo.plot.spectrum(line_bands=matched_bands)
# spec_malo.fit.frame(matched_bands, fit_conf=obs_cfg)
# spec_malo.plot.bands()

spec = hf2.get_spectrum(130, 117)
# spec.plot.spectrum()
spec.fit.frame(bands_file, cfgFile, line_detection=True)
spec.plot.spectrum()


# hf2.check.cube('H1_6563A', masks_file=spatial_mask_file)
# hf2.fit.spatial_mask(spatial_mask_file, bands=bands_file, fit_conf=obs_cfg, mask_list="MASK_0",
#                      line_detection=True, output_address=output_lines_log_file)

