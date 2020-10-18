import numpy as np

## ----------------  Rearrenge orientation -----------------------

matrix = np.array([[1, 2, 3, 8, 8], [0, 0, 4, 5, 6], [8, 8, 7, 8, 9], [1, 2, 0, 5, 2]])
X, Y = np.meshgrid(np.arange(matrix.shape[1]), np.arange(matrix.shape[0]))
X_flatten, Y_flatten = X.flatten(), Y.flatten()

print(matrix,'\n')
counter = 0
for idx_j in np.arange(matrix.shape[0]):
    for idx_i in np.arange(matrix.shape[1]):
        total_counter = np.ravel_multi_index((idx_j, idx_i), matrix.shape)
        quote = f'value {matrix[idx_j,idx_i]} == {matrix[Y_flatten[counter],X_flatten[counter]]} at == {matrix[Y_flatten[total_counter], X_flatten[total_counter]]}'
        quote2 = f' indeces: {idx_j}, {idx_i} == {np.unravel_index(counter, matrix.shape)}'
        print(quote, quote2)
        counter += 1

print(Y_flatten, '\n')
print(X_flatten, '\n')

print(matrix.size, Y_flatten.size, X_flatten.size)

## ----------------  Rearrenge orientation -----------------------

# matrix = np.array([[1,2,3,8,8],[0,0,4,5,6],[8,8,7,8,9],[1,2,0,5,2]])
# print('matrix original\n', matrix)
# matrix_1d = matrix.flatten()
# print('matrix flattend\n', matrix_1d)
# matrix_2d = B = np.reshape(matrix_1d, (4, 5))
# print('matrix reshaped\n', matrix_2d)



# ## ----------------  Dimension orientation -----------------------
#
# matrix = np.array([[1,2,3],[4,5,6],[7,8,9]])
# print(matrix)
# print(matrix[0, 2])
# print('First is y, second is x')

## ----------------  Squeeze -----------------------
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