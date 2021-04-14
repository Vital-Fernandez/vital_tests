import pandas as pd
from astropy.io import fits
from pathlib import Path

log_address = Path('/home/vital/Astro-data/Observations/gp121903_BR_linesLog_it1.txt')

linesDF = pd.read_csv(log_address, delim_whitespace=True, header=0, index_col=0)

sA = linesDF.to_records(index=True)

sA_types = sA.dtype.descr
for i, item in enumerate(sA_types):
    column, column_type = item
    if column_type == '|O':
        sA_types[i] = (column, 'U50')
sA = sA.astype(sA_types)

cols = fits.ColDefs(sA)
hdu=fits.BinTableHDU.from_columns(cols)
print(hdu.data['gauss_flux'])
print(hdu.data['index'][0])
print(type(hdu.data['index'][0]))


# linesDF.index = linesDF.index.map(str)
# sA = linesDF.loc[:, 'wavelength':'gauss_err'].to_records(index=True)
# index_Uni = sA['index'].astype('U25')
#
# sA_types = sA.dtype.descr
# sA_types[0] = ('index', 'U25')
# sA = sA.astype(sA_types)
#
# # for item in sA_types:
# #     column, column_type = item
# #     if column_type =
# # sA['index'] = index_Uni
#
# cols = fits.ColDefs(sA)
# hdu=fits.BinTableHDU.from_columns(cols)
# print(hdu.data['gauss_flux'])
# print(hdu.data['index'][0])
# print(type(hdu.data['index'][0]))