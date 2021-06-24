import numpy as np
import pandas as pd
from pathlib import Path

grid_columns = {'logZ': 'logZ',
                'logU': 'logU',
                'logNO': 'logNO',
                'carbon': 'carbon',
                'o3726': 'O2_3726A',
                'o3729': 'O2_3729A',
                'ne3869': 'Ne3_3869A',
                'ne3968': 'Ne3_3868A',
                'h3970': 'H1_3970A',
                's4070': 'S2_4069A',
                's4078': 'S2_4078A',
                'h4102': 'H1_4102A',
                'c4267': 'C2_4267A',
                'h4341': 'H1_4341A',
                'o4363': 'O3_4363A',
                'he4471': 'He1_4471A',
                'o4651': 'O2_4651A',
                'c4659': 'C2_4659A',
                'fe4668': 'Fe3_4668A',
                'he4686': 'He2_4686A',
                'ar4711': 'Ar4_4711A',
                'ar4740': 'Ar4_4740A',
                'h4861': 'H1_4861A',
                'o4959': 'O3_4959A',
                'o5007': 'O3_5007A',
                'ar5192': 'Ar3_5192A',
                'n5198': 'N1_5198A',
                'n5200': 'N1_5200A',
                'cl5518': 'Cl3_5518A',
                'cl5538': 'Cl3_5538A',
                'n5755': 'N2_5755A',
                'he5876': 'He1_5876A',
                'o6300': 'O1_6300A',
                's6312': 'S3_6312A',
                'n6548': 'N2_6548A',
                'h6563': 'H1_6563A',
                'n6584': 'N2_6584A',
                'he6678': 'He1_6678A',
                's6716': 'S2_6716A',
                's6731': 'S2_6731A',
                'he7065': 'He1_7065A',
                'ar7135': 'Ar3_7136A',
                'o7323': 'O2_7319A',
                'o7332': 'O2_7330A',
                'ar7751': 'Ar3_7751A',
                's9069': 'S3_9069A',
                's9532': 'H1_9229A'}

grid_file = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/HII-CHI-mistry_1Myr_grid.csv')
grid_DF = pd.read_csv(grid_file, skiprows=1, names=grid_columns.values())
grid_DF.logNO = np.round(grid_DF.logNO.values, decimals=3)


Halpha = np.power(10, grid_DF['H1_4341A'].values)

w = np.unique(grid_DF.logZ.values)
x = np.unique(grid_DF.logU.values)
y = np.unique(grid_DF.logNO.values)
z = np.unique(grid_DF.carbon.values)

file_matrix = grid_DF['O2_3726A'].values
file_mdimArray = file_matrix.reshape(len(w), len(x), len(y), len(z))

idx_w = 2
idx_x = 3
idx_y = 15
idx_z = 0
print(w[idx_w], x[idx_x], y[idx_y], z[idx_z])

print(file_mdimArray[idx_w, idx_x, idx_y, idx_z])
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