from pathlib import Path
import numpy as np
import pandas as pd
import lime

from astropy.io import fits
from astropy.table import Table
from desi_functions import desi_mask
from matplotlib import pyplot as plt

# Data for the new Sample files
output_columns = {'TARGETID':   'i8',
                  'SURVEY':     'U7',
                  'PROGRAM':    'U6',
                  'HEALPIX':    'i4',
                  'SPECTYPE':   'U6',
                  'SUBTYPE':    'U2',
                  'Z':          'f4',
                  'ZERR':       'f4',
                  'DELTACHI2':  'f4',
                  'ZCAT_PRIMARY': 'i1'}
input_columns = {k: 'category' if v.startswith('U') else v for k, v in output_columns.items()}

new_dict = {k: float(v) if isinstance(v, int) else v for k, v in original_dict.items()}


input_columns = {'TARGETID': 'Int64',
                 'SURVEY': 'category',
                 'PROGRAM': 'string',
                 'HEALPIX': 'Int32',
                 'SPECTYPE': 'string',
                 'SUBTYPE': 'string',
                 'Z': 'float32',
                 'ZERR': 'float32',
                 'DELTACHI2': 'float32',
                 'ZCAT_PRIMARY': 'Int8'}

column_names = list(output_columns.keys())

# Fits files with output targets
fits_sample_files = {'BGS_ANY':     'sample_BGS_ANY.fits',
                     'LRG':         'sample_LRG.fits',
                     'ELG':         'sample_ELG.fits',
                     'QSO':         'sample_QSO.fits',
                     'MWS_ANY':     'sample_MWS_ANY.fits',
                     'SCND_ANY':    'sample_SCND_ANY.fits'}

# Open DESI healpix master file
zpix_file = './zall-pix-fuji.fits'
fujidata = Table.read(zpix_file, hdu="ZCATALOG")

# Index DESI target
desi_target = fujidata["DESI_TARGET"]

check_dict = {}
for obj_type, output_file in fits_sample_files.items():

    # Get the masks which identify the object types
    obj_tgt_mask = desi_mask[obj_type]
    tgt_check = ((desi_target & obj_tgt_mask != 0) |
                 (fujidata["SV1_DESI_TARGET"] & obj_tgt_mask != 0) |
                 (fujidata["SV2_DESI_TARGET"] & obj_tgt_mask != 0) |
                 (fujidata["SV3_DESI_TARGET"] & obj_tgt_mask != 0))

    # Crop tables
    tb = fujidata[column_names][tgt_check]
    df = tb.to_pandas()

    assert len(tb) == df.index.size, f'Different number of rows'
    assert len(tb.columns) == df.columns.size, f'Different number of rows'

    # Save to fits
    lime.save_log(df, output_file, column_dtypes=output_columns)

    # Read the fits file:
    df2 = lime.load_log(output_file)

    # log_df = Table.read(output_file, 'LINELOG', character_as_bytes=False).to_pandas()
    # log_df.set_index('index', inplace=True)

    with fits.open(output_file) as hdul:
        hdu_data = hdul['LINELOG'].data
        df3 = pd.DataFrame(data=hdu_data)
        dtypes = [(key, value) for key, value in output_columns.items()]
        df4 = pd.DataFrame(data=hdu_data, dtype=dtypes)
        df4 = pd.DataFrame.from_records(hdu_data, index='index')
        df4 = df4.astype(input_columns)
        for column in df3:
            print(column.dtype)

df_memory = df4
memory_array_column = df_memory.memory_usage(deep=True)/(1024**2)
for i, column in enumerate(df_memory.columns):
    print(f'{column}) {df_memory[column].dtype} ({memory_array_column[column]:0.2f} MB)')


    # # Convert the types
    # for column in string_columns:
    #     df[column] = df[column].astype('category')
    #
    # for column in float32_columns:
    #     df[column] = df[column].astype('float32')

    # check_dict[obj_type] = tgt_check


# targets = list(fits_sample_files.keys())
# fig, ax = plt.subplots(1, 1, figsize=(8, 6))
#
# counts = []
# for tgt in targets:
#     n_targets = np.count_nonzero(check_dict[tgt])
#     counts.append(n_targets)
# ax.bar(targets, counts, color="purple", alpha=0.4, label="All spectra (includes duplicate targets)")
#
# for i in range(len(targets)):
#     ax.text(targets[i], counts[i], counts[i], ha="center", va="bottom", fontsize=16)
#
# ax.set_ylabel("Number of spectra")
# ax.semilogy()
# ax.set_ylim(1e5, 2e6)
#
# plt.legend(fontsize=18, frameon=False, markerfirst=False)
# plt.tight_layout()
# plt.show()


#
# (desi_target & elg_tgtmask != 0)
#
# for column in zpix_cat.columns:
#     print(column, zpix_cat[column].dtype)
#
# target_columns = ['TARGETID', 'SURVEY', 'PROGRAM', 'HEALPIX', 'Z', 'ZERR', 'DELTACHI2', 'SPECTYPE', 'SUBTYPE']
# string_columns = ['SURVEY', 'PROGRAM', 'SPECTYPE', 'SUBTYPE']
# float32_columns = ['Z', 'ZERR', 'CHI2']
#
# tb = zpix_cat[target_columns]
# df = zpix_cat[target_columns].to_pandas()
#
# assert len(tb) == df.index.size, f'Different number of rows'
# assert len(tb.columns) == df.columns.size, f'Different number of rows'
#
#
# memory_array_column = df.memory_usage(deep=True)/(1024**2)
# for i, column in enumerate(df.columns):
#     print(f'{column}) {df[column].dtype} ({memory_array_column[column]:0.2f} MB)')
#
#
# for column in string_columns:
#     df[column] = df[column].astype('category')
#
# for column in float32_columns:
#     df[column] = df[column].astype('float32')
#
#
# memory_array_column2 = df.memory_usage(deep=True)/(1024**2)
# for i, column in enumerate(df.columns):
#     print(f'{column}) {df[column].dtype} ({memory_array_column[column]:0.2f} -> {memory_array_column2[column]:0.2f} MB)')
#
# df2 = df.set_index(['TARGETID', 'SURVEY', 'PROGRAM'])
# memory_array_column3 = df2.memory_usage(deep=True)/(1024**2)
#
# params_dtype = {'TARGETID': 'i8', 'SURVEY': 'U7', 'PROGRAM': 'U6',
#                 'HEALPIX': 'i4', 'Z': 'f4', 'ZERR': 'f4', 'DELTACHI2': 'f4',
#                 'SPECTYPE': 'U6', 'SUBTYPE': 'U2'}
#
# lime.save_log(df, 'single_index.fits', column_dtypes=params_dtype)
# # lime.save_log(df2, 'multi_index.fits')
#
# df3 = lime.load_log('single_index.fits')
# tb2 = Table.read('single_index.fits', 'LINELOG', character_as_bytes=False)
# df3 = tb2.to_pandas()
# df3.set_index('index', inplace=True)
#
#
# tb3 = Table.read('single_index.fits', 'LINELOG', character_as_bytes=True)
# df4 = tb3.to_pandas()

# params_dtype = {'TARGETID': np.dtype('b'), 'SURVEY': str, 'PROGRAM': str,
#                 'HEALPIX': np.int32, 'Z': np.float32, 'ZERR': np.float32, 'CHI2': np.float32,
#                 'SPECTYPE': str, 'SUBTYPE': str}
#
# linesSA = log.to_records(index=True, column_dtypes=params_dtype, index_dtypes=np.int64)