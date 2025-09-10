import theano
import numpy as np


# def myOperation(a_value, b_value):
#     return a_value**2 + 2 * a_value * b_value
#
#
# a = theano.tensor.vector()
# b = theano.tensor.dscalar()
# myOperation_tt = theano.function([a, b], myOperation)
# # print(myOperation_tt(a_value=np.array([1, 2, 3]), b_value=2))

a = theano.tensor.vector() # declare variable
b = theano.tensor.dscalar()
out = a ** 2 + 2 * a * b              # build symbolic expression
f = theano.function([a, b], out)   # compile function
print(f([0, 1, 2], 3))