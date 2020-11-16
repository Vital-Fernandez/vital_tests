import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
import matplotlib.pyplot as plt

conf_file_address = 'gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)
dataFolder = Path(obsData['file_information']['data_folder'])

files_dict = dict(spec_mio = dataFolder/'gp004054_BR.fits',
                    spec_ric = dataFolder/'GP004054_comb_new_clean.fits',
                    spec_ric_resf = dataFolder/'GP004054_comb_new_rf_clean.fits')

z = 0.283232

wave1, flux_array1, header1 = sr.import_fits_data(files_dict['spec_mio'], instrument='OSIRIS')
wave2, flux_array2, header2 = sr.import_fits_data(files_dict['spec_ric'], instrument='OSIRIS')
wave3, flux_array3, header3 = sr.import_fits_data(files_dict['spec_ric_resf'], instrument='OSIRIS')

fig, ax = plt.subplots()
ax.step(wave2 * 1/(1+z), flux_array2, label='spec_mio')
ax.step(wave3, flux_array3, label='spec_ric_resf', linestyle=':')
ax.legend()
plt.show()
ric_idx_max = np.argmax(flux_array2)
mio_idx_max = np.argmax(flux_array3)
print(flux_array2[ric_idx_max])
print(flux_array3[mio_idx_max])

print(flux_array1/flux_array2)
print((flux_array1-flux_array2).sum())


