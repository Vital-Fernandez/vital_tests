import numpy as np
import pandas as pd
import inspect
from pathlib import Path
from matplotlib import pyplot as plt, rcParams
from lime.plots import STANDARD_PLOT

rcParams.update(STANDARD_PLOT)


def get_values(df, idcs, column):
    return df.loc[idcs, column].values


def slice_grid(param_x, param_y, grid_df, **kwargs):

    # Container for the slice indeces
    idcs_slice = np.full(grid_df.index.size, True)

    # Loop through the grid parameter axes
    for ax_param, ax_value in kwargs.items():
        if ax_param != param_x:
            if ax_param in grid_df.columns:
                idcs_slice = idcs_slice & (grid_df[ax_param] == ax_value)
            else:
                print(f'-WARNING: {ax_param} parameter not found in grid')

    x_value = grid_df.loc[idcs_slice, param_x].values
    y_values = grid_df.loc[idcs_slice, param_y].values

    return x_value, y_values


grid_columns_large = {'logZ': 'logOH',
                    'logU': 'logU',
                    'logNO': 'logNO',
                    'carbon': 'carbon',
                    'o3726': 'O2_3726A',
                    'o3729': 'O2_3729A',
                    'ne3869': 'Ne3_3869A',
                    'ne3968': 'Ne3_3968A',
                    'h3970': 'H1_3970A',
                    's4070': 'S2_4069A',
                    's4078': 'S2_4078A',
                    'h4102': 'H1_4102A',
                    'c4267': 'C2_4267A',
                    'h4341': 'H1_4341A',
                    'o4363': 'O3_4363A',
                    'he4471': 'He1_4471A',
                    'o4651': 'O1_4651A',
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
                    's9532': 'S3_9531A'}

grid_columns_HIICHI_mistry = {'12+log(O/H)':'logOH',
                             'log(N/O)':'logNO',
                             'log(U)':'logU',
                             'OII_3727':'O2_3726A_b',
                             'NeIII_3868':'Ne3_3869A',
                             'OIII_4363':'O3_4363A',
                             'OIII_5007':'O3_5007A',
                             'NII_6584':'N2_6584A',
                             'SII_6717,31':'S2_6716A_b'}

line_labels = {'O3_5007A': r'Flux $[OIII]5007\AA$',
               'N2_6584A': r'Flux $[NII]6584\AA$',
               'O2_7319A': r'Flux $[OII]7319\AA$',
               'O2_7319A_b': r'Flux $[OII]7319,7330\AA$',
               'O2_3726A_b': r'Flux $[OII]3726,3729\AA$',
               'logOH': '12+log(O/H)',
               'logNO': 'log(N/O)',
               'logU': 'log(U)'}

# grid_files_dict = {'grid_Jorge': Path('D:/Dropbox/Astrophysics/Papers/muse_CGCG007/data/HII-CHI-mistry_1Myr_grid.csv'),
#                    'grid_Epm': Path('D:/Dropbox/Astrophysics/Papers/muse_CGCG007/data/formated_log_C17_Popstar_1Myr.dat'),
#                    'HIIchimistry_v2': Path('D:/Dropbox/Astrophysics/Tools/HCm_v2.0/C13_cha_1Myr_v2.0.dat'),
#                    'HIIchimistry_v3': Path('D:/Dropbox/Astrophysics/Tools/HCm_v3.0/C17_cha_1Myr_v3.1.dat'),
#                    'HIIchimistry_v4': Path('D:/Dropbox/Astrophysics/Tools/HCm_v4.2/C17_cha_1Myr_v4.0.dat'),
#                    'HIIchimistry_v5': Path('D:/Dropbox/Astrophysics/Tools/HCm_v5.22/Libraries_opt/C17_POPSTAR_1myr.dat')}
grid_files_dict = {'grid_Jorge': Path('/home/vital/Dropbox/Astrophysics/Papers/muse_CGCG007/data/HII-CHI-mistry_1Myr_grid.csv'),
                   'grid_Epm': Path('/home/vital/Dropbox/Astrophysics/Papers/muse_CGCG007/data/formated_log_C17_Popstar_1Myr.dat'),
                   'HIIchimistry_v2': Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v2.0/C13_cha_1Myr_v2.0.dat'),
                   'HIIchimistry_v3': Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v3.0/C17_cha_1Myr_v3.1.dat'),
                   'HIIchimistry_v4': Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v4.2/C17_cha_1Myr_v4.0.dat'),
                   'HIIchimistry_v5': Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v5.22/Libraries_opt/C17_POPSTAR_1myr.dat')}

grid_column_dict = {'grid_Jorge': grid_columns_large,
                    'grid_Epm': {},
                    'HIIchimistry_v2': grid_columns_HIICHI_mistry,
                    'HIIchimistry_v3': grid_columns_HIICHI_mistry,
                    'HIIchimistry_v4': grid_columns_HIICHI_mistry,
                    'HIIchimistry_v5': grid_columns_HIICHI_mistry}


# Load the data
grid_jorge = pd.read_csv(grid_files_dict['grid_Jorge'], skiprows=1, names=grid_columns_large.values())
grid_jorge.iloc[:, 4:] = np.power(10, grid_jorge.iloc[:, 4:])
grid_jorge['O2_3726A_b'] = grid_jorge['O2_3726A'] + grid_jorge['O2_3729A']
grid_jorge['O2_7319A_b'] = grid_jorge['O2_7319A'] + grid_jorge['O2_7330A']

grid_epm = pd.read_csv(grid_files_dict['grid_Epm'], delim_whitespace=True)
grid_epm.iloc[:, 3:] = np.power(10, grid_epm.iloc[:, 3:])
grid_epm['O2_3726A_b'] = grid_epm['O2_3726A'] + grid_epm['O2_3729A']
grid_epm['O2_7319A_b'] = grid_epm['O2_7319A'] + grid_epm['O2_7330A']

file_columns = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'OII_3727', 'OIII_4363', 'OIII_5007', 'NII_6584', 'SII_6717,31']
HIIchimistry_v2 = pd.read_csv(grid_files_dict['HIIchimistry_v2'], delim_whitespace=True, names=file_columns)
HIIchimistry_v2.rename(columns=grid_columns_HIICHI_mistry, inplace=True)

file_columns = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'OII_3727', 'OIII_4363', 'OIII_5007', 'NII_6584', 'SII_6717,31']
HIIchimistry_v3 = pd.read_csv(grid_files_dict['HIIchimistry_v3'], delim_whitespace=True, names=file_columns)
HIIchimistry_v3.rename(columns=grid_columns_HIICHI_mistry, inplace=True)

file_columns = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'OII_3727', 'NeIII_3868', 'OIII_4363', 'OIII_5007', 'NII_6584', 'SII_6717,31']
HIIchimistry_v4 = pd.read_csv(grid_files_dict['HIIchimistry_v4'], delim_whitespace=True, names=file_columns, comment='#')
HIIchimistry_v4.rename(columns=grid_columns_HIICHI_mistry, inplace=True)

HIIchimistry_v5 = pd.read_csv(grid_files_dict['HIIchimistry_v5'], delim_whitespace=True)
HIIchimistry_v5.rename(columns=grid_columns_HIICHI_mistry, inplace=True)

# Container for the grids
grid_dict = dict(grid_jorge=grid_jorge,
                 grid_epm=grid_epm,
                 HIIchimistry_v2=HIIchimistry_v2,
                 HIIchimistry_v3=HIIchimistry_v3,
                 HIIchimistry_v4=HIIchimistry_v4,
                 HIIchimistry_v5=HIIchimistry_v5,)

# Testing points
logOH_ref, logU_ref, logNO_ref = 7.8, -2.5, -1.5

# OIII fluxes
param_list = ['logOH', 'logNO', 'logU', 'logNO', 'logU']
line_list = ['O3_5007A', 'N2_6584A', 'N2_6584A', 'O2_3726A_b', 'O2_3726A_b']

# # Plot the intervals
# for param, line in zip(param_list, line_list):
#
#     fig, ax = plt.subplots()
#     x_jorge_O, y_jorge_O = slice_grid(param, line, grid_jorge, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='O')
#     x_jorge_N, y_jorge_N = slice_grid(param, line, grid_jorge, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='N')
#     x_epm, y_epm = slice_grid(param, line, grid_epm, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
#
#     x_HIIchim_v2, y_HIIchim_v2 = slice_grid(param, line, HIIchimistry_v2, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
#     x_HIIchim_v3, y_HIIchim_v3 = slice_grid(param, line, HIIchimistry_v3, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
#     x_HIIchim_v4, y_HIIchim_v4 = slice_grid(param, line, HIIchimistry_v4, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
#     x_HIIchim_v5, y_HIIchim_v5 = slice_grid(param, line, HIIchimistry_v5, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
#
#     ax.scatter(x_jorge_O, y_jorge_O, label='Jorge Grid carbon = O', marker='*')
#     ax.scatter(x_jorge_N, y_jorge_N, label='Jorge Grid carbon = N', marker='s')
#     ax.scatter(x_epm, y_epm, label='Enrique Grid', marker='^')
#
#     ax.scatter(x_HIIchim_v2, y_HIIchim_v2, label='HII-CHI-mistry v2.00 grid', alpha=0.5)
#     # ax.scatter(x_HIIchim_v3, y_HIIchim_v3, label='HII-CHI-mistry v3.00 grid', alpha=0.5)
#     # ax.scatter(x_HIIchim_v4, y_HIIchim_v4, label='HII-CHI-mistry v4.20 grid', alpha=0.5)
#     ax.scatter(x_HIIchim_v5, y_HIIchim_v5, label='HII-CHI-mistry v5.22 grid', alpha=0.5)
#
#     ax.legend()
#     ax.update({'xlabel': f'{line_labels[param]}', 'ylabel': line_labels[line]})
#     plt.show()

# # NII fluxes versus N/O for various O
# param_list = ['logNO', 'logNO']
# line_list = ['N2_6584A', 'O2_3726A_b']
#
# for param, line in zip(param_list, line_list):
#
#     fig, ax = plt.subplots()
#
#     for OH in [7.6, 7.8, 8.0]:
#         x_jorge_O, y_jorge_O = slice_grid(param, line, grid_jorge, logOH=OH, logU=logU_ref, logNO=logNO_ref, carbon='O')
#         x_jorge_N, y_jorge_N = slice_grid(param, line, grid_jorge, logOH=OH, logU=logU_ref, logNO=logNO_ref, carbon='N')
#         x_epm, y_epm = slice_grid(param, line, grid_epm, logOH=OH, logU=logU_ref, logNO=logNO_ref)
#
#         x_HIIchim_v2, y_HIIchim_v2 = slice_grid(param, line, HIIchimistry_v2, logOH=OH, logU=logU_ref, logNO=logNO_ref)
#         x_HIIchim_v3, y_HIIchim_v3 = slice_grid(param, line, HIIchimistry_v3, logOH=OH, logU=logU_ref, logNO=logNO_ref)
#         x_HIIchim_v4, y_HIIchim_v4 = slice_grid(param, line, HIIchimistry_v4, logOH=OH, logU=logU_ref, logNO=logNO_ref)
#         x_HIIchim_v5, y_HIIchim_v5 = slice_grid(param, line, HIIchimistry_v5, logOH=OH, logU=logU_ref, logNO=logNO_ref)
#
#         # ax.scatter(x_csv, y_csv, label='CSV grid Carbon = O', marker='*')
#         # ax.scatter(x_csv_N, y_csv_N, label='CSV grid Carbon = N', marker='*')
#         # ax.scatter(x_HIIchim_v2, y_HIIchim_v2, label='HII-CHI-mistry v2.00 grid', alpha=0.5)
#         # ax.scatter(x_HIIchim_v3, y_HIIchim_v3, label='HII-CHI-mistry v3.00 grid', alpha=0.5)
#         # ax.scatter(x_HIIchim_v4, y_HIIchim_v4, label='HII-CHI-mistry v4.20 grid', alpha=0.5)
#         ax.scatter(x_HIIchim_v5, y_HIIchim_v5, label=f'12 + log(O/H) = {OH} (jorge)', alpha=0.5)
#         ax.scatter(x_epm, y_epm, label=f'12 + log(O/H) = {OH} (epm)', alpha=0.5,  marker='*')
#
#     ax.legend()
#     ax.update({'xlabel': f'{line_labels[param]}', 'ylabel': line_labels[line]})
#     plt.show()

# NII fluxes versus N/O for various O
param_list = ['logNO', 'logNO']
line_list = ['N2_6584A', 'O2_3726A_b']

for param, line in zip(param_list, line_list):

    fig, ax = plt.subplots()

    for OH in [7.6, 7.8, 8.0]:


        x_epm, y_epm = slice_grid(param, line, grid_epm, logOH=OH, logU=logU_ref, logNO=logNO_ref)

        idcs_below = x_epm <= -1.25
        print(f'Flux deviation at O/H = {OH} is {y_epm[idcs_below].std()}')
        print(len(y_epm[idcs_below]), np.sum(idcs_below))
        print()
        x_HIIchim_v5, y_HIIchim_v5 = slice_grid(param, line, HIIchimistry_v5, logOH=OH, logU=logU_ref, logNO=logNO_ref)
        ax.scatter(x_epm, y_epm, label=f'12 + log(O/H) = {OH} (epm)', alpha=0.5,  marker='*')

    ax.legend()
    ax.update({'xlabel': f'{line_labels[param]}', 'ylabel': line_labels[line]})
    plt.show()

