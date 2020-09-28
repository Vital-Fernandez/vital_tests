import numpy as np

## ----------------  Dimension orientation -----------------------

matrix = np.array([[1,2,3],[4,5,6],[7,8,9]])
print(matrix)
print(matrix[0, 2])
print('First is y, second is x')

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