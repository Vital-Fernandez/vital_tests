from pathlib import Path
import numpy as np
import pandas as pd
import lime

from astropy.io import fits
from astropy.table import Table


zpix_file = './zall-pix-fuji.fits'
zpix_cat = Table.read(zpix_file, hdu="ZCATALOG")
surveys = np.unique(zpix_cat['SURVEY'])
program = np.unique(zpix_cat['PROGRAM'])
spectype = np.unique(zpix_cat['SPECTYPE'])
subtype = np.unique(zpix_cat['SUBTYPE'])


for column in zpix_cat.columns:
    print(column, zpix_cat[column].dtype)

target_columns = ['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'Z', 'ZERR', 'DELTACHI2', 'SPECTYPE', 'SUBTYPE']
string_columns = ['SURVEY', 'PROGRAM', 'SPECTYPE', 'SUBTYPE']
float32_columns = ['Z', 'ZERR', 'CHI2']

tb = zpix_cat[target_columns]
df = zpix_cat[target_columns].to_pandas()

assert len(tb) == df.index.size, f'Different number of rows'
assert len(tb.columns) == df.columns.size, f'Different number of rows'


memory_array_column = df.memory_usage(deep=True)/(1024**2)
for i, column in enumerate(df.columns):
    print(f'{column}) {df[column].dtype} ({memory_array_column[column]:0.2f} MB)')


for column in string_columns:
    df[column] = df[column].astype('category')

for column in float32_columns:
    df[column] = df[column].astype('float32')


memory_array_column2 = df.memory_usage(deep=True)/(1024**2)
for i, column in enumerate(df.columns):
    print(f'{column}) {df[column].dtype} ({memory_array_column[column]:0.2f} -> {memory_array_column2[column]:0.2f} MB)')

df2 = df.set_index(['TARGETID', 'SURVEY', 'PROGRAM'])
memory_array_column3 = df2.memory_usage(deep=True)/(1024**2)

params_dtype = {'TARGETID': 'i8', 'SURVEY': 'U7', 'PROGRAM': 'U6',
                'HEALPIX': 'i4', 'Z': 'f4', 'ZERR': 'f4', 'DELTACHI2': 'f4',
                'SPECTYPE': 'U6', 'SUBTYPE': 'U2'}

lime.save_log(df, 'single_index.fits', column_dtypes=params_dtype)
# lime.save_log(df2, 'multi_index.fits')

df3 = lime.load_log('single_index.fits')
tb2 = Table.read('single_index.fits', 'LINELOG', character_as_bytes=False)
df3 = tb2.to_pandas()
df3.set_index('index', inplace=True)


tb3 = Table.read('single_index.fits', 'LINELOG', character_as_bytes=True)
df4 = tb3.to_pandas()

# params_dtype = {'TARGETID': np.dtype('b'), 'SURVEY': str, 'PROGRAM': str,
#                 'HEALPIX': np.int32, 'Z': np.float32, 'ZERR': np.float32, 'CHI2': np.float32,
#                 'SPECTYPE': str, 'SUBTYPE': str}
#
# linesSA = log.to_records(index=True, column_dtypes=params_dtype, index_dtypes=np.int64)