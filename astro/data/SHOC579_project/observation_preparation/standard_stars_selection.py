from pathlib import Path
import numpy as np
import src.specsiser as sr
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
from matplotlib import pyplot as plt

conf_file = Path('../obsConf.ini')
obsData = sr.loadConfData(conf_file)

data_folder = Path(obsData['data_location']['data_folder'])
results_folder = Path(obsData['data_location']['results_folder'])

osirisStdStars_file = '/home/vital/Dropbox/Astrophysics/TelescopeTime/Standard_Stars/osirisStandardStars.csv'

esoStdStars_file = '/home/vital/Dropbox/Astrophysics/TelescopeTime/Standard_Stars/esoStandardStars.csv'

osirisStdStarsDF = pd.read_csv(osirisStdStars_file, header=0, index_col=0)
esoStdStarsDF = pd.read_csv(esoStdStars_file, header=0, index_col=0)

star_list = osirisStdStarsDF.index.values
ra_array, dec_array = np.zeros(len(star_list)), np.zeros(len(star_list))
for i, star in enumerate(star_list):
    ra = coord.Angle(osirisStdStarsDF.loc[star, 'RA'], unit=u.hourangle)
    dec = coord.Angle(osirisStdStarsDF.loc[star, 'DEC'], unit=u.degree)
    ra_array[i], dec_array[i] = ra.radian, dec.radian
    print(i, star, ra_array[i], ra, dec_array[i], dec, dec.degree)

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="aitoff")
ax.scatter(ra_array, dec_array)
ax.grid(True)
ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
plt.show()


