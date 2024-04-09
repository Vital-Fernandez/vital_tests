import lime
import numpy as np

spec_address = './Hf2_2.fits'
mask_address = './Hf2_2_mask.fits'
mask_plot = './Hf2_2_mask.png'
hf2 = lime.Cube.from_file(spec_address, instrument='MUSE')

pctls = np.array([75])
hf2.spatial_masking('H1_6563A', param='SN_line', contour_pctls=pctls, output_address=mask_address)
hf2.plot.cube('H1_6563A', masks_file=mask_address, output_address=mask_plot)


