import src.specsiser as sr
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# data_folder = Path('D:/Dropbox/Astrophysics/Data/J083843/blue')
# file_name = 'combined_2d_apall.fits'

data_folder = Path('D:/Dropbox/')
file_name = 'J0838_blue_apall.fits'


wave, flux_array, header = sr.import_fits_data(data_folder/file_name, instrument='ISIS', frame_idx=0)

print(header)

labelsDict = {'xlabel': r'Wavelength $(\AA)$',
              'ylabel': r'Flux $(erg\,cm^{-2} s^{-1} \AA^{-1})$',
              'title': f'Galaxy J083843 extracted blue arm'}

print(np.sum(flux_array[0][0]))

# Plot spectra components
fig, ax = plt.subplots(figsize=(12, 8))
ax.plot(wave, flux_array[2][0], label='Object flux')
ax.update(labelsDict)
ax.legend()
plt.show()
