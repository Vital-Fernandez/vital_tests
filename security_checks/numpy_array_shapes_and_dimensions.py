import numpy as np


def function_a(a, b, c, d):
    return a + b + c + d/100000


w = np.arange(6.0, 9.0, 0.5)
x = np.arange(1, 1000, 200)
y = np.arange(-4.0, -1.5, 0.5)
z = np.arange(1e6, 1e7, 1e6)

file_matrix = np.zeros((w.size*x.size*y.size*z.size, 5))
i = 0
for w_value in w:
    for x_value in x:
        for y_value in y:
            for z_value in z:
                file_matrix[i, 0] = w_value
                file_matrix[i, 1] = x_value
                file_matrix[i, 2] = y_value
                file_matrix[i, 3] = z_value
                file_matrix[i, 4] = function_a(w_value, x_value, y_value, z_value)
                i += 1

file_mdimArray = file_matrix[:, 4].reshape(len(w), len(x), len(y), len(z))

print(file_mdimArray[0, 0, 0, 0])
print(function_a(w[0], x[0], y[0], z[0]))

# w_grid, x_grid, y_grid, z_grid = np.meshgrid(w, x, y, z)
#
# fa = function_a(w_grid, x_grid, y_grid, z_grid)
# fb = function_b(w_grid, x_grid, y_grid, z_grid)
#
# print(fa)

# a = np.array([[1, 2, 3, 4, 5],
#               [5, 6, 7, 8, 9]])
#
# print(a.mean(axis=1))
# print(a[0, :].mean(axis=1))

# print(a)
# print(a.shape)
# print(a.ndim)
# print()
# print(a[:, 1])
# print(a[:, 1].shape)
# print(a[:, 1].ndim)
# print()
# print(a[0, 1])
# print(a[0, 1].shape)
# print(a[0, 1].ndim)
# print()



# def myFunc(a_in):
#     a_in = np.array(a_in, ndmin=1)
#     return a_in
#
# a = 5
# print('Antes', a)
# b = myFunc(a)
# print('Despues', a, b)

## ----------------  Ndimensions -----------------------
# a = np.array([[6645.603569, 6663.100335, 6711.874023, 6735.247559, 6750.335676, 6766.603503],
#              [7027.765924, 7049.636882, 7062.26123, 7075.280273, 7089.317048, 7102.127181]])
# b = np.array([6645.603569, 6663.100335, 6711.874023, 6735.247559, 6750.335676, 6766.603503])
# c = np.array([[6645.603569, 6663.100335, 6711.874023, 6735.247559, 6750.335676, 6766.603503]])
#
# a1 = np.array(a, ndmin=2)
# b1 = np.array(b, ndmin=2)
# c1 = np.array(c, ndmin=2)
#
# a2 = np.squeeze(a1)
# b2 = np.squeeze(b1)
# c2 = np.squeeze(c1)
#
# print(a.shape)
# print(b.shape)
# print(c.shape)
# print()
# print(a1.shape)
# print(b1.shape)
# print(c1.shape)
# print()
# print(a2.shape)
# print(b2.shape)
# print(c2.shape)


## ----------------  Rearrenge orientation -----------------------

# matrix = np.array([[1, 2, 3, 8, 8], [0, 0, 4, 5, 6], [8, 8, 7, 8, 9], [1, 2, 0, 5, 2]])
# X, Y = np.meshgrid(np.arange(matrix.shape[1]), np.arange(matrix.shape[0]))
# X_flatten, Y_flatten = X.flatten(), Y.flatten()
#
# print(matrix,'\n')
# counter = 0
# for idx_j in np.arange(matrix.shape[0]):
#     for idx_i in np.arange(matrix.shape[1]):
#         total_counter = np.ravel_multi_index((idx_j, idx_i), matrix.shape)
#         quote = f'value {matrix[idx_j,idx_i]} == {matrix[Y_flatten[counter],X_flatten[counter]]} at == {matrix[Y_flatten[total_counter], X_flatten[total_counter]]}'
#         quote2 = f' indeces: {idx_j}, {idx_i} == {np.unravel_index(counter, matrix.shape)}'
#         print(quote, quote2)
#         counter += 1
#
# print(Y_flatten, '\n')
# print(X_flatten, '\n')
#
# print(matrix.size, Y_flatten.size, X_flatten.size)

## ----------------  Rearrenge orientation -----------------------

# matrix = np.array([[1,2,3,8,8],[0,0,4,5,6],[8,8,7,8,9],[1,2,0,5,2]])
# print('matrix original\n', matrix)
# matrix_1d = matrix.flatten()
# print('matrix flattend\n', matrix_1d)
# matrix_2d = np.reshape(matrix_1d, (4, 5))
# print('matrix reshaped\n', matrix_2d)



# ## ----------------  Dimension orientation -----------------------
#
# matrix = np.array([[1,2,3],[4,5,6],[7,8,9]])
# print(matrix)
# print(matrix[0, 2])
# print('First is y, second is x')

# # ----------------  Squeeze -----------------------
# a = np.array([1, 2, 3, 4, 5])
# b = np.array(5)
# c = np.array(5, ndmin=1)
# d = np.array([1, 2, 3, 4, 5], ndmin=1)
# e = np.array(d, ndmin=1)
# f = np.array([7])
#
# print(a, a.shape)
# print(b, b.shape, type(b))
# print(c, c.shape)
# print(d, d.shape)
# print(e, e.shape)
# print(f, f.shape)
#
# coso = 'H1_4863A'
# caso = np.array(coso.split('-'), ndmin=1)
# print(caso, caso.shape)
#
# coso = 'H1_4863A-H1_6563A'
# caso = np.array(coso.split('-'), ndmin=1)
# print(caso, caso.shape)
#
# print('\n')
# a = np.empty(3)
# b = np.empty((3, 1))
# print(a)
# print(b)
# print(b.squeeze())
# print(b)
# b = b.squeeze()
# print(b)