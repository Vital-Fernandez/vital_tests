from pathlib import Path
import numpy as np
import src.specsiser as sr
import pandas as pd
import astropy.coordinates as coord
import astropy.units as u
from matplotlib import pyplot as plt
import mplcursors


megara_library_stars = ['HR1544', 'HR5501', 'HR5501', 'HR5501', 'HR1544']

conf_file = Path('../obsConf.ini')
obsData = sr.loadConfData(conf_file)

data_folder = Path(obsData['data_location']['data_folder'])
results_folder = Path(obsData['data_location']['results_folder'])

targets_file = '/home/vital/Dropbox/Astrophysics/TelescopeTime/MEGARA_NOexcess_Gradients_2021/phase2/targets_coordinates.txt'
osirisStdStars_file = '/home/vital/Dropbox/Astrophysics/TelescopeTime/Standard_Stars/osirisStandardStars.csv'
esoStdStars_file = '/home/vital/Dropbox/Astrophysics/TelescopeTime/Standard_Stars/esoStandardStars.csv'

osirisStdStarsDF = pd.read_csv(osirisStdStars_file, header=0, index_col=0)
osirisStdStarsDF['type'] = 'Osiris Stars'

esoStdStarsDF = pd.read_csv(esoStdStars_file, header=0, index_col=0)
esoStdStarsDF['type'] = 'ESO Stars'

targetsDF = pd.read_csv(targets_file, delim_whitespace=True, header=0, index_col=0)
targetsDF['type'] = 'Targets'

objDF = pd.concat([osirisStdStarsDF, esoStdStarsDF, targetsDF])
objDF['RA_hour'] = coord.Angle(objDF['RA'], unit=u.hourangle)
objDF['DEC_degree'] = coord.Angle(objDF['DEC'], unit=u.degree)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
objList = []
for i, type in enumerate(np.unique(objDF['type'].values)):
    idcs = (objDF.type == type)
    labels = objDF.loc[idcs].index.values
    ra_array, dec_array = objDF.loc[idcs, 'RA_hour'].values, objDF.loc[idcs, 'DEC_degree'].values
    data = ax.scatter(ra_array, dec_array, label=type)
    mplcursors.cursor(data).connect("add", lambda sel, labels=labels: sel.annotation.set_text(labels[sel.target.index]))

ax.set_xlim(0, 24)
ax.set_ylim(-90, 90)
ax.update({'xlabel': r'RA (hr)', 'ylabel': r'Dec (deg)', 'title': 'Standard Stars vs Objects'})
ax.legend()
plt.show()


MEGARA_file = '/home/vital/Dropbox/Astrophysics/TelescopeTime/MEGARA_NOexcess_Gradients_2021/phase2/Appendix_C_MEGASTAR_MEGARA_GTC_Stellar_Spectral_Library_II_First_Release.csv'
targets_file = '/home/vital/Dropbox/Astrophysics/TelescopeTime/MEGARA_NOexcess_Gradients_2021/phase2/targets_coordinates.txt'

MegaraStarsDF = pd.read_csv(MEGARA_file, delimiter=';', header=0, index_col=0)
MegaraStarsDF['type'] = 'MEGARA Stars'

targetsDF = pd.read_csv(targets_file, delim_whitespace=True, header=0, index_col=0)
targetsDF['type'] = 'Targets'

objDF = pd.concat([MegaraStarsDF, targetsDF])

objDF['RA_hour'] = coord.Angle(objDF['RA'], unit=u.hourangle)
objDF['DEC_degree'] = coord.Angle(objDF['DEC'], unit=u.degree)

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111)
objList = []
for type in np.unique(objDF.type.values):
    idcs = (objDF.type == type)
    labels = objDF.loc[idcs].index.values
    ra_array, dec_array = objDF.loc[idcs, 'RA_hour'].values, objDF.loc[idcs, 'DEC_degree'].values
    data = ax.scatter(ra_array, dec_array, label=type)
    mplcursors.cursor(data).connect("add", lambda sel, labels=labels: sel.annotation.set_text(labels[sel.target.index]))

ax.set_xlim(0, 24)
ax.set_ylim(-90, 90)
ax.update({'xlabel': r'RA (hr)', 'ylabel': r'Dec (deg)', 'title': 'Standard Stars vs Objects'})
ax.legend()
plt.show()