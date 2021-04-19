import pandas as pd
from astropy.io import fits
from pathlib import Path
import numpy as np
from astropy.table import Table

log_address = Path('/home/vital/Astro-data/Observations/gp121903_BR_linesLog_it1.txt')
fits_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/treatment/CGCG007/CGCG007_linesLog.fits')

# default_types = {'index': '<U50'),
#                  ('ion': '<U50'),
#                  ('latexLabel', '<U50'),
#                  ('blended', '<U50'),
#                  ('observation', '<U50'),
#                  ('comments', '<U50')])

default_types = {'index': '<U50',
                 'wavelength': '<f8',
                 'intg_flux': '<f8',
                 'intg_err': '<f8',
                 'gauss_flux': '<f8',
                 'gauss_err': '<f8',
                 'eqw': '<f8',
                 'eqw_err': '<f8',
                 'ion': '<U50',
                 'pynebCode': '<f8',
                 'pynebLabel': '<f8',
                 'lineType': '<f8',
                 'latexLabel': '<U50',
                 'blended': '<U50',
                 'w1': '<f8',
                 'w2': '<f8',
                 'w3': '<f8',
                 'w4': '<f8',
                 'w5': '<f8',
                 'w6': '<f8',
                 'm_continuum': '<f8',
                 'n_continuum': '<f8',
                 'cont': '<f8',
                 'std_continuum': '<f8',
                 'peak_flux': '<f8',
                 'peak_wave': '<f8',
                 'snr_line': '<f8',
                 'snr_cont': '<f8',
                 'amp': '<f8',
                 'mu': '<f8',
                 'sigma': '<f8',
                 'amp_err': '<f8',
                 'mu_err': '<f8',
                 'sigma_err': '<f8',
                 'v_r': '<f8',
                 'v_r_err': '<f8',
                 'sigma_vel': '<f8',
                 'sigma_err_vel': '<f8',
                 'observation': '<U50',
                 'comments': '<U50',
                 'obsFlux': '<f8',
                 'obsFluxErr': '<f8',
                 'f_lambda': '<f8',
                 'obsInt': '<f8',
                 'obsIntErr': '<f8'}

linesDF = pd.read_csv(log_address, delim_whitespace=True, header=0, index_col=0)
linesSA = linesDF.to_records(index=True, column_dtypes=default_types, index_dtypes='<U50')
linesCol = fits.ColDefs(linesSA)
linesHDU = fits.BinTableHDU.from_columns(linesCol, name='lineLog')
linesFits = Path('/home/vital/Astro-data/Observations/test_HDU.fits')
linesHDU.writeto(linesFits, overwrite=True)

linesData = fits.getdata(linesFits, 'lineLog', ver=1)
linesDF_reload = Table.read(linesFits, 'lineLog', character_as_bytes=False).to_pandas()
linesDF_reload.set_index('index', inplace=True)

# linesDF = pd.read_csv(log_address, delim_whitespace=True, header=0, index_col=0)
#
# sA = linesDF.to_records(index=True)
#
# sA_dict = dict(sA.dtype.descr)
# sA_dict.update(default_types)
# sA = sA.astype(list(sA_dict.items()))
#
#
# # sA_dict = dict(zip(sA.dtype.names, sA.dtype.descr))
# # for i, item in enumerate(sA_types):
# #     column, column_type = item
# #     if column_type == '|O':
# #         sA_types[i] = (column, 'U50')
# # sA = sA.astype(sA_types)
# #
# #
# # lst2.extend(item[0] for item in lst)
# #
# cols = fits.ColDefs(sA)
# extname = f'151-151_linelog'
# hdu = fits.BinTableHDU.from_columns(cols, name=extname)
#
# print(hdu)

# fits.update(fits_address, data=hdu.data, header=hdu.header, extname=extname)
# fits.append(fits_address, data=hdu.data, header=hdu.header, extname=extname)
# fits_data = fits.info(fits_address)
# hdr1 = fits.getheader(fits_address, extname=extname)
# hdul_list = fits.open(fits_address)

# extname = '170-170_LINELOG'
# linesBinHDU = fits.getdata(fits_address, extname, ver=1)
# linesSA = np.array(linesBinHDU)
# linesDF = pd.DataFrame(linesSA)
#
# linesTable = Table.read(fits_address, extname)
# linesDF = Table.read(fits_address, extname).to_pandas()
# linesDF.set_index('index', inplace=True)
# linesTable.to_pandas()
# print(linesDF)




# from astropy.utils.data import get_pkg_data_filename
# from astropy.table import Table
# from astropy.io import fits
#
# fits_file = get_pkg_data_filename('tutorials/FITS-tables/chandra_events.fits')
#
# fits.info(fits_file)





