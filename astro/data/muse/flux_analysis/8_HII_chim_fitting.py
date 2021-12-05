import numpy as np
import pandas as pd
from pathlib import Path
import time
from astropy.io import fits
from astropy.table import Table
from astro.data.muse.common_methods import grid_HII_CHI_mistry_conversion as labelConver
from astro.papers.gtc_greenpeas.common_methods import epm_HII_CHI_mistry


# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

for i, obj in enumerate(objList):

    if i == 0:

        # Data location
        cube_address = fitsFolder/fileList[i]
        objFolder = resultsFolder/obj
        voxelFolder = resultsFolder/obj/'voxel_data'
        db_address = objFolder / f'{obj}_database.fits'
        fitsLog_address = objFolder / f'{obj}_linesLog.fits'

        # Output data
        HIIchimistry_fits = objFolder / f'{obj}_HIIchimistry.fits'
        hdul_lineslog = fits.HDUList()

        # Loop throught the line regions
        start = time.time()
        for idx_region in [0, 1, 2]:

            region_label = f'region_{idx_region}'
            region_mask = fits.getdata(db_address, region_label, ver=1)
            region_mask = region_mask.astype(bool)
            idcs_voxels = np.argwhere(region_mask)

            # Loop through the region voxels
            n_voxels = idcs_voxels.shape[0]
            for idx_voxel, idx_pair in enumerate(idcs_voxels):

                print(f'-- Treating voxel {idx_voxel}/{n_voxels} ({idx_pair})')
                idx_j, idx_i = idx_pair
                logLabel = f'{idx_j}-{idx_i}_linelog'
                # output_fit_file = objFolder/'voxel_treatments'/f'{idx_j}-{idx_i}_HII_CHI_mistry_fit.txt'

                # Load voxel lines log
                linesLog_BinTable = fits.getdata(fitsLog_address, logLabel, ver=1)
                linesDF = Table(linesLog_BinTable).to_pandas()
                linesDF.set_index('index', inplace=True)

                # Prepare data for HII-CHI-mistry
                idcs_inputLines = linesDF.index.isin(labelConver.keys())
                input_lines = linesDF.loc[idcs_inputLines].index.values

                HII_CHI_mistry_DF = pd.DataFrame()
                HII_CHI_mistry_DF.loc[0, 'ID'] = logLabel
                flux_Hbeta = linesDF.loc['H1_4861A', 'intg_flux']
                for lineLabel in input_lines:
                    HII_CHI_mistry_DF.loc[0, labelConver[lineLabel]] = linesDF.loc[lineLabel, 'intg_flux'] / flux_Hbeta
                    HII_CHI_mistry_DF.loc[0, f'e{labelConver[lineLabel]}'] = linesDF.loc[lineLabel, 'intg_err'] / flux_Hbeta
                lineSA = HII_CHI_mistry_DF.to_records(index=False) #column_dtypes=default_linelog_types, index_dtypes='<U50')

                # Run HII-CHI-mistry
                outputSA = epm_HII_CHI_mistry(lineSA, output_file='None', n=200, sed=1, inter=1)
                linesCol = fits.ColDefs(outputSA)
                linesHDU = fits.BinTableHDU.from_columns(linesCol, name=f'{idx_j}-{idx_i}_HIIchimistry')
                hdul_lineslog.append(linesHDU)

            # Store the drive
            hdul_lineslog.writeto(HIIchimistry_fits, overwrite=True, output_verify='fix')
        end = time.time()

print(f'- Execution time {(end - start)/60:.3f} min')

# grid_file = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/HII-CHI-mistry_1Myr_grid.csv')
# grid_DF = pd.read_csv(grid_file, skiprows=1, names=grid_columns.values())
# grid_DF.logNO = np.round(grid_DF.logNO.values, decimals=3)
#
# Halpha = np.power(10, grid_DF['H1_4341A'].values)
#
# w = np.unique(grid_DF.logZ.values)
# x = np.unique(grid_DF.logU.values)
# y = np.unique(grid_DF.logNO.values)
# z = np.unique(grid_DF.carbon.values)
#
# file_matrix = grid_DF['O2_3726A'].values
# file_mdimArray = file_matrix.reshape(len(w), len(x), len(y), len(z))
#
# idx_w = 2
# idx_x = 3
# idx_y = 15
# idx_z = 0
# print(w[idx_w], x[idx_x], y[idx_y], z[idx_z])
#
# print(file_mdimArray[idx_w, idx_x, idx_y, idx_z])
# import pandas as pd
# import pymysql
# import matplotlib.pyplot as plt
# import os
#
# os.environ['MdB_HOST'] = '3mdb.astro.unam.mx'
# os.environ['MdB_USER'] = 'oiii5007'
# os.environ['MdB_PASSWD'] = 'OVN_user'
# os.environ['MdB_PORT'] = '3306'
# os.environ['MdB_DB_17'] = '3MdB_17'
# os.environ['MdB_DBs'] = '3MdBs'
# os.environ['MdB_DBp'] = '3MdB'
#
# # co = pymysql.connect(host=os.environ['MdB_HOST'],
# #                      db=os.environ['MdB_DB_17'],
# #                      user=os.environ['MdB_USER'],
# #                      passwd=os.environ['MdB_PASSWD'])
# co = pymysql.connect(host='3mdb.astro.unam.mx', db='3MdB_17', user='OVN_user', passwd='oiii5007')
# # co = pymysql.connect(host='3mdb.astro.unam.mx',
# #                      db='3MdB_17',
# #                      user='oiii5007',
# #                      passwd='OVN_user')
# res = pd.read_sql("select log10(N__2_658345A/H__1_656281A) as n2, log10(O__3_500684A/H__1_486133A) as o3, OXYGEN as O from tab_17 where ref = 'BOND'", con=co)
# co.close()
#
# f, ax = plt.subplots()
# sc = ax.scatter(res['n2'], res['o3'], c=12+res['O'], edgecolor='')
# ax.set_xlabel("log [NII]/Ha")
# ax.set_ylabel("log [OIII]/Hb")
# cb = f.colorbar(sc, ax=ax)
# cb.set_label("12 + log O/H")
# plt.show()