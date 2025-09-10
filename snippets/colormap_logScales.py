from pylab import *
import numpy as np
from matplotlib.colors import SymLogNorm
from matplotlib import pyplot as plt

np.random.seed(10)

data = np.random.uniform(low=-10, high=10, size=(10,10))
print(np.mean(data))
norm = SymLogNorm(5, vmin=-10, vmax=10)
fig, axes = plt.subplots()
im = axes.imshow(data,extent=[-10,10,-10,10],cmap=plt.cm.jet,norm=norm)
cb = fig.colorbar(im)
plt.show()