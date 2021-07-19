from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import src.specsiser as sr

c_KMpS = 299792.458

obsData = sr.loadConfData('./xshooter_LzLCS.ini')
data_folder = Path(obsData['data_location']['data_folder'])
results_folder = Path(obsData['data_location']['results_folder'])
objfile_list = obsData['data_location']['objfile_list']
sigmafile_list = obsData['data_location']['sigmafile_list']
objRef_list = obsData['data_location']['ref_list']
maskfile = obsData['data_location']['generalMask']

wmin_array = obsData['sample_data']['w_min_array']
wmax_array = obsData['sample_data']['w_max_array']
norm_flux = obsData['sample_data']['norm_flux']
z_obj = obsData['sample_data']['z_obj']
profile_conf = obsData['line_fitting']

DF_list = [None, None]
for i, objName in enumerate(objRef_list):

    # input data
    lineslog_file = results_folder/f'{objName}_linesLog.txt'

    # Load data
    linesDF = sr.lineslogFile_to_DF(lineslog_file)
    DF_list[i] = linesDF

# Join the dataframes
objDF = pd.concat(DF_list)

print(objDF)

theo_waves = objDF.wavelength.values
obs_waves = objDF.center.values

objDF['z_gauss'] = obs_waves/theo_waves - 1

idcs_wide = objDF.index.str.contains('_w1')
idcs_narrow = ~idcs_wide & ~objDF.index.str.contains('_w2')

objDF.loc[idcs_narrow, 'z_gauss']
objDF.loc[idcs_wide, 'z_gauss']

z_narrow = objDF.loc[idcs_narrow, 'z_gauss'].values
z_wide = objDF.loc[idcs_wide, 'z_gauss'].values
labels = [r'Narrow $z_{{median}} = {{{:0.6f}}}$'.format(np.median(z_narrow)),
          r'Wide $z_{{median}} = {{{:0.6f}}}$'.format(np.median(z_wide))]

fig, ax = plt.subplots()
a = ax.boxplot([z_narrow, z_wide], labels=labels)
ax.scatter(np.ones(z_narrow.size), z_narrow, alpha=0.2)
ax.scatter(np.ones(z_wide.size)*2, z_wide, alpha=0.2)
ax.update({'xlabel':'Line component', 'ylabel':'Line Redshift', 'title':'Redshift for wide and narrow components'})
plt.tight_layout()
plt.show()

v_narrow = objDF.loc[idcs_narrow, 'z_gauss'].values * c_KMpS
v_wide = objDF.loc[idcs_wide, 'z_gauss'].values * c_KMpS
labels = [r'Narrow $v_{{median}} = {{{:0.0f}}}\,km/s$'.format(np.median(v_narrow)),
          r'Wide $v_{{median}} = {{{:0.0f}}}\,km/s$'.format(np.median(v_wide))]

fig, ax = plt.subplots()
ax.boxplot([v_narrow, v_wide], labels=labels)
ax.scatter(np.ones(v_narrow.size), v_narrow, alpha=0.2)
ax.scatter(np.ones(v_wide.size)*2, v_wide, alpha=0.2)
ax.update({'xlabel':'Line component', 'ylabel':'Line Redshift (km/s)', 'title':'Redshift for wide and narrow components'})
plt.tight_layout()
plt.show()

output_log = results_folder/f'line_measurements_z_gaussian.txt'
output_columns = ['wavelength', 'gauss_flux', 'gauss_err', 'eqw', 'z_gauss']
with open(output_log, 'w') as f:
    f.write(objDF.loc[:, output_columns].to_string(index=True, index_names=False))
