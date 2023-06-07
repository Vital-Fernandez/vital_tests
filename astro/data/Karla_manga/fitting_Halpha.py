import numpy as np

import numpy as np
import lime
import pandas as pd
from pathlib import Path


path_fits = Path('NGC5471_C_34-28-LR-R_sky.txt')
configFile = Path('ejemplo_ajuste.cfg')

obs_cfg = lime.load_cfg(configFile)

logN2 = lime.load_log(f'LR-V_log.txt')
logHbeta = lime.load_log(f'LR-B_log.txt')
log_df = result = pd.concat([logHbeta, logN2])

wave_LR_R, flux_LR_R = np.loadtxt(path_fits, unpack=True, usecols=(0, 1))
band_edges = np.array([6438.03, 6508.66, 6525.10, 6616.0, 6627.70, 6661.82])

gp_spec = lime.Spectrum(wave_LR_R, flux_LR_R, redshift=0.00091, norm_flux=0.00001)
gp_spec.load_log(log_df)
gp_spec.plot.spectrum(rest_frame=True)

line = 'H1_6563A_b'
gp_spec.fit.band(line, band_edges, fit_conf=obs_cfg['LR-R_line_fitting'])
gp_spec.plot.band()


#
#
#
# # gp_spec.plot.spectrum(label='NGC5471', rest_frame=True)
#
# # line = 'H1_6562.7192A_b'
# # gp_spec.fit.band(line, band_edges, fit_conf=obs_cfg['manga_2compsHR_plus_wide_line_fitting'])
# # gp_spec.plot.band()
# # gp_spec.save_log(f'manga_2compsHR_plus_wide_line_fitting.txt')
#
# # line = 'H1_6562.7192A_b'
# # gp_spec.fit.band(line, band_edges, fit_conf=obs_cfg['manga_3compsHR_line_fitting'])
# # gp_spec.plot.band()
# #
# line = 'H1_6562.7192A_b'
# gp_spec.fit.band(line, band_edges, fit_conf=obs_cfg['auroral_import_line_fitting'])
# gp_spec.plot.band()
# gp_spec.save_log('LR-R_log.txt')
#
# print(lime.__version__)