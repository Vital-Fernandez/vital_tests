import numpy as np
import lime

from astropy.io import fits
lime.show_instrument_cfg()

fname = '/home/vital/Astrodata/LzLCS_MIRI/release/data/s3d_cubes/J003601+003307_ch2.fits'
obj = lime.Cube.from_file(fname, instrument='miri', redshift=0.35)
# obj.check.cube('H1_74599A', rest_frame=True)
lime.show_instrument_cfg()

fname = '/home/vital/Astrodata/LzLCS_MIRI/release/data/x1d_pipe_channel/J003601+003307_pipe1_ch2_x1d.fits'
obj = lime.Spectrum.from_file(fname, instrument='nirspec', redshift=0.35)
# obj.plot.spectrum()

lime.show_instrument_cfg()


# TTYPE1  = 'WAVELENGTH'
# TFORM1  = 'D       '
# TUNIT1  = 'um      '
# TTYPE2  = 'FLUX    '
# TFORM2  = 'D       '
# TUNIT2  = 'Jy      '
# TTYPE3  = 'FLUX_ERROR'
# TFORM3  = 'D       '
# TUNIT3  = 'Jy      '