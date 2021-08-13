import sys
import numpy as np
import matplotlib.pyplot as plt

x = np.array([1, 2, 3, 4])
y = np.array([1, 2, 3, 4])

label = 'python {}.{}'.format(sys.version_info.major, sys.version_info.minor)

fig, ax = plt.subplots()
ax.plot(x, y, label=label)
ax.legend()
ax.update({'xlabel':'Wavelength', 'ylabel':'Flux', 'title':'Gaussian fitting'})
plt.show()
