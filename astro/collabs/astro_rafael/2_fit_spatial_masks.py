import numpy as np
import lime

# spec_address = './Hf2_2.fits'
# conf_file = './muse.toml'
# mask_address = './Hf2_2_mask_2.fits'
# bands_file = './Lineas_muse_minimo.txt'
# measurements_file = './Hf2_2_line_measurements.fits'
# cube = lime.Cube.from_file(spec_address, instrument='MUSE', redshift=0)
# cube.fit.spatial_mask(mask_address, bands=bands_file, fit_conf=conf_file,
#                       line_detection=True, output_address=measurements_file)
# cube.check.cube('H1_6563A', lines_file=measurements_file, masks_file=mask_address)

spec_address = './NGC3242.fits'
conf_file = 'muse.toml'
mask_address = './NGC3242_mask.fits'
bands_file = './Lineas_muse_minimo.txt'
measurements_file = './NGC3242_line_measurements.fits'
cube = lime.Cube.from_file(spec_address, instrument='MUSE', redshift=0)
cube.fit.spatial_mask(mask_address, bands=bands_file, fit_conf=conf_file,
                      line_detection=True, output_address=measurements_file)
cube.check.cube('H1_6563A', lines_file=measurements_file, masks_file=mask_address)


# spec_address = './M1-42_30S.fits'
# conf_file = './muse.toml'
# mask_address = './M1-42_30S_mask.fits'
# bands_file = './Lineas_muse_minimo.txt'
# measurements_file = './M1-42_30S_line_measurements.fits'
# cube = lime.Cube.from_file(spec_address, instrument='MUSE', redshift=0)
#
# cube.fit.spatial_mask(mask_address, bands=bands_file, fit_conf=conf_file,
#                       line_detection=True, output_address=measurements_file)
# cube.check.cube('H1_6563A', lines_file=measurements_file, masks_file=mask_address)
