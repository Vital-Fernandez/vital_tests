import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model

def Funk(x, k, y0, good_points=None):  # note: add keyword argument
    f = k*x + y0
    if good_points is not None:
        f = f[good_points]       # apply mask of good data points
    return f

x = np.array([1,2,3,4,5, 6,7,8.,9,10])
y = np.array([1,2,3,4,5,30,35.,40,45,50])
y += np.random.normal(size=len(x), scale=0.19) # add some noise to make it fun

# make an array of the indices of the "good data points"
# does not need to be contiguous.
good_points=np.array([0,1,2,3,4])

# turn your model function Funk into an lmfit Model
mymodel = Model(Funk)

# create parameters, giving initial values. Note that parameters are
# named using the names of your function's argument and that keyword
# arguments with non-numeric defaults like 'good points' are seen to
#  *not* be parameters. Like the independent variable `x`, you'll
# need to pass that in when you do the fit.
# also: parameters can be fixed, or given `min` and `max` attributes

params = mymodel.make_params(k=1.4,  y0=0.2)
params['k'].min = 0

# do the fit to the 'good data', passing in the parameters, the
# independent variable `x` and the `good_points` mask.
print(type(mymodel))
result = mymodel.fit(y[good_points], params, x=x, good_points=good_points)

# print out a report of best fit values, uncertainties, correlations, etc.
print(result.fit_report())

# plot the results, again using the good_points array as needed.
plt.plot(x, y, 'o', label='all data')
plt.plot(x[good_points], result.best_fit[good_points], label='fit to good data')
plt.legend()
plt.show()