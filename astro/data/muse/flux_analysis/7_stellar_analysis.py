# from pathlib import Path
# from astropy.io import fits
#
# # file_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/CGCG007/Papaderos_20210616/red_cgcg007025_lowZChab_HBIN024_FDres2.fits')
# # file_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/CGCG007/Papaderos_20210616/red_cgcg007025_3D_SLres_HBIN024.fits')
# # file_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/CGCG007/Papaderos_20210616/red_cgcg007025_HBIN024_FDres2.fits')
# file_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/CGCG007/Papaderos_20210616/red_u5205_HBIN030_FDres2.fits')
#
# hdul = fits.open(file_address)
# print(fits.info(file_address))
# hdul.close()

import pyneb as pn
from matplotlib import pyplot as plt

# Notebooks tutoriales PyNeb https://github.com/Morisset/PyNeb_devel/tree/master/docs

O2 = pn.Atom('O', 2)

diags = pn.Diagnostics()
OIII_4363_5007_ratio = 0.02
OII_3726_3729_ratio = 1.00

TOIII, neOII = diags.getCrossTemDen('[OIII] 4363/5007', '[OII] 3726/3729',
                                    OIII_4363_5007_ratio, OII_3726_3729_ratio)
print(TOIII, neOII)

O2_EG = pn.EmisGrid('O', 2, n_tem=30, n_den=30)
f, ax = plt.subplots(figsize=(7, 5))
O2_EG.plotContours(to_eval='L(3726)/L(3729)', ax=ax, log_levels=False)
plt.show()
