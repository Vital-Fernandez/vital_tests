import numpy as np

# obsWave = np.array([4981.62, 5030.54, 6593.03, 9111.79])
# trueWaves = np.array([4959, 5007, 6563, 9069])
trueWaves = np.array([4341.471, 4861, 4958.911, 5006.843])
obsWave = np.array([5093.6, 5705, 5819.64, 5876.13])

z_array = obsWave/trueWaves

z_mean, z_std = z_array.mean(), z_array.std()

print(z_array)
print(z_mean, z_std)
print(obsWave / (z_mean))
