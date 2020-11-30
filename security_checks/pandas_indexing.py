import numpy as np
import pandas as pd

num_array = np.array([[1, 1, 3, 0, 2],
                      [4, 2, 6, 0, 8],
                      [4, 3, 9, 5, 7],
                      [5, 1, -5, 3, np.nan]])

sample_df = pd.DataFrame(data=num_array, columns=['a', 'b', 'c', 'd', 'e'])

# Object database as a pandas dataframe
db_headers = ['y_voxel', 'x_voxel', 'voxel_label', 'flux']
size_y, size_x = num_array.shape[0], num_array.shape[1]
obsDF = pd.DataFrame(index=np.arange(size_y * size_x), columns=db_headers)

# Voxel indeces where Cube shape = (lambda , Y, X)
y_range, x_range = np.ogrid[0:num_array.shape[0], 0:num_array.shape[1]]
X, Y = np.meshgrid(x_range, y_range)
X_flatten, Y_flatten = X.flatten(), Y.flatten()
obsDF['y_voxel'] = Y_flatten
obsDF['x_voxel'] = X_flatten
obsDF['flux'] = num_array.flatten()

voxel_label_array = np.empty(obsDF.index.size, dtype=str)
voxel_label_array.fill('-')
voxel_label_array = np.core.defchararray.add(X_flatten.astype(str), voxel_label_array)
voxel_label_array = np.core.defchararray.add(voxel_label_array, Y_flatten.astype(str))
obsDF['voxel_label'] = voxel_label_array

idx_pair = (1, 2)
true_value = num_array[idx_pair]

idx_j, idx_i = idx_pair
idx_db = (obsDF.y_voxel == idx_j) & (obsDF.x_voxel == idx_i)

print(true_value, obsDF.loc[idx_db, 'flux'].values)


print(obsDF)



# # Voxel indeces where Cube shape = (lambda , Y, X)
# y_range, x_range = np.ogrid[0:num_array.shape[0], 0:num_array.shape[1]]
# Y, X = np.meshgrid(x_range, y_range)
# X_flatten, Y_flatten = X.flatten(), Y.flatten()
# obsDF['y_voxel'] = X_flatten
# obsDF['x_voxel'] = Y_flatten
# obsDF['flux'] = num_array.flatten()
#
# y_range, x_range = np.ogrid[0:num_array.shape[0], 0:num_array.shape[1]]
# Y, X = np.meshgrid(y_range, x_range)
# voxel_label_array = np.empty(obsDF.index.size, dtype=str)
# voxel_label_array.fill('-')
# voxel_label_array = np.core.defchararray.add(X_flatten.astype(str), voxel_label_array)
# voxel_label_array = np.core.defchararray.add(voxel_label_array, Y_flatten.astype(str))
# obsDF['voxel_label'] = voxel_label_array
#
# print(obsDF)