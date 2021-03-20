import numpy as np

Teff_array = np.array([57620, 56591, 56894])
TeffErr_array = np.array([2960, 3178, 2116])

print(Teff_array.mean(), TeffErr_array.mean())

Teff_array = np.array([62129, 60600, 61566])
TeffErr_array = np.array([4387, 2824, 2018])

print(Teff_array.mean(), TeffErr_array.mean())


Teff_array = np.array([74335, 65070, 70342])
TeffErr_array = np.array([536, 334, 238])

print(Teff_array.mean(), TeffErr_array.mean())

logU_array = np.array([2.55, 2.16, 2.13])
logUErr_array = np.array([0.01, 0.01, 0.01])

print(logU_array.mean(), logUErr_array.mean())
