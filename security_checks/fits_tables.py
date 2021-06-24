import pandas as pd
from astropy.io import fits
from pathlib import Path
import numpy as np
from astropy.table import Table
import src.specsiser as sr

# log_address = Path('/home/vital/Astro-data/Observations/gp121903_BR_linesLog_it1.txt')
# fits_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/treatment/CGCG007/CGCG007_linesLog.fits')
#
# default_types = {'index': '<U50',
#                  'wavelength': '<f8',
#                  'intg_flux': '<f8',
#                  'intg_err': '<f8',
#                  'gauss_flux': '<f8',
#                  'gauss_err': '<f8',
#                  'eqw': '<f8',
#                  'eqw_err': '<f8',
#                  'ion': '<U50',
#                  'pynebCode': '<f8',
#                  'pynebLabel': '<f8',
#                  'lineType': '<f8',
#                  'latexLabel': '<U50',
#                  'blended_label': '<U50',
#                  'w1': '<f8',
#                  'w2': '<f8',
#                  'w3': '<f8',
#                  'w4': '<f8',
#                  'w5': '<f8',
#                  'w6': '<f8',
#                  'm_cont': '<f8',
#                  'n_cont': '<f8',
#                  'cont': '<f8',
#                  'std_cont': '<f8',
#                  'peak_flux': '<f8',
#                  'peak_wave': '<f8',
#                  'snr_line': '<f8',
#                  'snr_cont': '<f8',
#                  'amp': '<f8',
#                  'mu': '<f8',
#                  'sigma': '<f8',
#                  'amp_err': '<f8',
#                  'mu_err': '<f8',
#                  'sigma_err': '<f8',
#                  'v_r': '<f8',
#                  'v_r_err': '<f8',
#                  'sigma_vel': '<f8',
#                  'sigma_err_vel': '<f8',
#                  'observation': '<U50',
#                  'comments': '<U50',
#                  'obsFlux': '<f8',
#                  'obsFluxErr': '<f8',
#                  'f_lambda': '<f8',
#                  'obsInt': '<f8',
#                  'obsIntErr': '<f8'}
#
# linesDF = pd.read_csv(log_address, delim_whitespace=True, header=0, index_col=0)
# linesSA = linesDF.to_records(index=True, column_dtypes=default_types, index_dtypes='<U50')
# linesCol = fits.ColDefs(linesSA)
# linesHDU = fits.BinTableHDU.from_columns(linesCol, name='lineLog')
# linesFits = Path('/home/vital/Astro-data/Observations/test_HDU.fits')
# linesHDU.writeto(linesFits, overwrite=True)
#
# linesData = fits.getdata(linesFits, 'lineLog', ver=1)
# linesDF_reload = Table.read(linesFits, 'lineLog', character_as_bytes=False).to_pandas()
# linesDF_reload.set_index('index', inplace=True)

outputDb = '/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/treatment/CGCG007/voxel_data/166-166_fitting.db'
db_fits = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/treatment/CGCG007/CGCG007_db.fits')

fit_results = sr.load_MC_fitting(outputDb)

total_params_list = np.array(list(fit_results['Fitting_results'].keys()))
inputLabels = fit_results['Input_data']['lineLabels_list']
inputFlux = fit_results['Input_data']['inputFlux_array']
inputErr = fit_results['Input_data']['inputErr_array']

# Save variable traces
list_columns = []
for i, param in enumerate(fit_results['trace'].varnames):
    trace = fit_results['trace'][param]
    if param in total_params_list:
        col_param = fits.Column(name=param, format='E', array=trace)
        list_columns.append(col_param)
    else:
        if ('_Op' not in param) and ('_log__' not in param):
            col_param = fits.Column(name=param, format='E', array=trace)
            list_columns.append(col_param)

# Save flux traces
if 'calcFluxes_Op' in fit_results['trace'].varnames:
    trace = fit_results['trace']['calcFluxes_Op']
    for i in range(trace.shape[1]):
        flux_trace = trace[:, i]
        col_param = fits.Column(name=inputLabels[i], format='E', array=flux_trace)
        list_columns.append(col_param)

# Generate BinTable
extname = f'{151}-{153}_linelog'
cols = fits.ColDefs(list_columns)
hdu = fits.BinTableHDU.from_columns(cols, name=extname)

# Safe extra data in header
for param, value in fit_results['Fitting_results'].items():
    hdu.header[f'hierarch {param}'] = value[0]
    hdu.header[f'hierarch {param}_err'] = value[1]

for i, label in enumerate(inputLabels):
    hdu.header[f'hierarch flux_{label}'] = inputFlux[i]
    hdu.header[f'hierarch err_{label}'] = inputErr[i]

fits.update(db_fits, data=hdu.data, header=hdu.header, extname=extname, verify=True)

if db_fits.is_file():
    try:
        print('Updating')
        fits.update(db_fits, data=hdu.data, header=hdu.header, extname=extname, verify=True)
    except:
        print('Appending')
        fits.append(db_fits, data=hdu.data, header=hdu.header, extname=extname)
else:
    print('Writting')
    hdu.writeto(db_fits, overwrite=True, output_verify='fix')

print(fits.info(db_fits))
print(fits.getheader(db_fits, extname))
data = fits.getdata(db_fits, extname)
hdr = fits.getheader(db_fits, extname)
print(hdr['flux_H1_6563A'])
print(data['cHbeta'])
print(data['T_low'])

