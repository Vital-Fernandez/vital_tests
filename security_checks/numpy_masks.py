import numpy as np
from astropy.io import fits

A = np.matrix([[1, 2, 3, 33],
               [4, 5, 6, 66],
               [7, 8, 9, 99]])

# Los valores de la mascara (true) es donde no se realizan operaciones con ese array
A = np.ma.masked_where(A > 6, A)
A_mask = A.mask.astype(int)

print('\nBoolean mask')
print(A.mask)

print('\nInteger mask')
print(A_mask)

idcs_voxels = np.argwhere(A_mask)
n_voxels = idcs_voxels.shape[0]

print('\nOriginal mask')
for idx_voxel in np.arange(n_voxels):
    idx_j, idx_i = idcs_voxels[idx_voxel]
    print(idx_j, idx_i, A_mask.data[idx_j, idx_i])

# Saving the file
hdul_masks = fits.HDUList()
hdul_masks.append(fits.PrimaryHDU())
mask_hdu = fits.ImageHDU(name='Mask_extenstion', data=A_mask, ver=1)
hdul_masks.append(mask_hdu)
hdul_masks.writeto('mask_example.fits', overwrite=True, output_verify='fix')

# Loading the file
region_mask = fits.getdata('mask_example.fits', 'Mask_extenstion')
idcs_voxels = np.argwhere(region_mask)
n_voxels = idcs_voxels.shape[0]

print('\nLoaded mask')
for idx_voxel in np.arange(n_voxels):
    idx_j, idx_i = idcs_voxels[idx_voxel]
    print(idx_j, idx_i, A_mask.data[idx_j, idx_i])

