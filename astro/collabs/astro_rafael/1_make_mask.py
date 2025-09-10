import lime
import numpy as np


spec_address = './Hf2_2.fits'
mask_address = './Hf2_2_mask_2.fits'
bands_file = './Lineas_muse_minimo.txt'
cube = lime.Cube.from_file(spec_address, instrument='MUSE', redshift=0)
pctls = np.array([75])
cube.spatial_masking('H1_6563A', bands=bands_file, param='SN_line', contour_pctls=pctls, output_address=mask_address)
cube.check.cube('H1_6563A', bands=bands_file, rest_frame=True, masks_file=mask_address)

spec_address = './NGC3242.fits'
mask_address = './NGC3242_mask.fits'
bands_file = './Lineas_muse_minimo.txt'
pctls = np.array([53])
cube = lime.Cube.from_file(spec_address, instrument='MUSE', redshift=0)
cube.spatial_masking('H1_6563A', bands=bands_file, param='SN_line', contour_pctls=pctls, output_address=mask_address)
cube.check.cube('H1_6563A', bands=bands_file, rest_frame=True, masks_file=mask_address)

spec_address = './M1-42_30S.fits'
mask_address = './M1-42_30S_mask.fits'
bands_file = './Lineas_muse_minimo.txt'
pctls = np.array([75])
cube = lime.Cube.from_file(spec_address, instrument='MUSE', redshift=0)
cube.spatial_masking('H1_6563A', bands=bands_file, param='SN_line', contour_pctls=pctls, output_address=mask_address)
cube.check.cube('H1_6563A', bands=bands_file, rest_frame=True, masks_file=mask_address)
