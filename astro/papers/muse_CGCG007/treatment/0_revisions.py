import numpy as np
import pandas as pd
import inspect
from pathlib import Path
from matplotlib import pyplot as plt, rcParams
from lime.plots import STANDARD_PLOT

rcParams.update(STANDARD_PLOT)


def get_values(df, idcs, column):
    return df.loc[idcs, column].values


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

# Load the data
large_grid_file = Path('/home/vital/Dropbox/Astrophysics/Papers/muse_CGCG007/data/HII-CHI-mistry_1Myr_grid.csv')
large_DF = pd.read_csv(large_grid_file, skiprows=1, names=grid_columns_large.values())

large_DF.iloc[:, 4:] = np.power(10, large_DF.iloc[:, 4:])
large_DF['O2_7319A_b'] = large_DF['O2_7319A'] + large_DF['O2_7330A']
large_DF['O2_3726A_b'] = large_DF['O2_3726A'] + large_DF['O2_3729A']

HIICHImistry_file = Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v2.0/C13_cha_1Myr_v2.0.dat')
file_columns = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'OII_3727', 'OIII_4363', 'OIII_5007', 'NII_6584', 'SII_6717,31']
HIICHIm_DF_v2 = pd.read_csv(HIICHImistry_file, delim_whitespace=True, names=file_columns)
HIICHIm_DF_v2.rename(columns=grid_columns_HIICHI_mistry, inplace=True)

HIICHImistry_file = Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v3.0/C17_cha_1Myr_v3.1.dat')
file_columns = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'OII_3727', 'OIII_4363', 'OIII_5007', 'NII_6584', 'SII_6717,31']
HIICHIm_DF_v3 = pd.read_csv(HIICHImistry_file, delim_whitespace=True, names=file_columns)
HIICHIm_DF_v3.rename(columns=grid_columns_HIICHI_mistry, inplace=True)

HIICHImistry_file = Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v4.2/C17_cha_1Myr_v4.0.dat')
file_columns = ['12+log(O/H)', 'log(N/O)', 'log(U)', 'OII_3727', 'NeIII_3868', 'OIII_4363', 'OIII_5007', 'NII_6584', 'SII_6717,31']
HIICHIm_DF_v4 = pd.read_csv(HIICHImistry_file, delim_whitespace=True, names=file_columns, comment='#')
HIICHIm_DF_v4.rename(columns=grid_columns_HIICHI_mistry, inplace=True)

HIICHImistry_file = Path('/home/vital/Dropbox/Astrophysics/Tools/HCm_v5.22/Libraries_opt/C17_POPSTAR_1myr.dat')
HIICHIm_DF_v5 = pd.read_csv(HIICHImistry_file, delim_whitespace=True)
HIICHIm_DF_v5.rename(columns=grid_columns_HIICHI_mistry, inplace=True)

# Testing points
logOH_ref, logU_ref, logNO_ref = 7.8, -2.5, -1.5

# OIII fluxes
param_list = ['logOH', 'logNO', 'logU', 'logNO', 'logU']
line_list = ['O3_5007A', 'N2_6584A', 'N2_6584A', 'O2_3726A_b', 'O2_3726A_b']

# Plot the intervals
for param, line in zip(param_list, line_list):
    fig, ax = plt.subplots()
    x_csv, y_csv = slice_grid(param, line, large_DF, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='O')
    x_csv_N, y_csv_N = slice_grid(param, line, large_DF, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='N')
    x_HIIchim_v5, y_HIIchim_v5 = slice_grid(param, line, HIICHIm_DF_v5, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
    x_HIIchim_v2, y_HIIchim_v2 = slice_grid(param, line, HIICHIm_DF_v2, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
    x_HIIchim_v3, y_HIIchim_v3 = slice_grid(param, line, HIICHIm_DF_v3, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
    x_HIIchim_v4, y_HIIchim_v4 = slice_grid(param, line, HIICHIm_DF_v4, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)

    ax.scatter(x_csv, y_csv, label='CSV grid Carbon = O', marker='*')
    ax.scatter(x_csv_N, y_csv_N, label='CSV grid Carbon = N', marker='*')
    ax.scatter(x_HIIchim_v2, y_HIIchim_v2, label='HII-CHI-mistry v2.00 grid', alpha=0.5)
    ax.scatter(x_HIIchim_v3, y_HIIchim_v3, label='HII-CHI-mistry v3.00 grid', alpha=0.5)
    ax.scatter(x_HIIchim_v4, y_HIIchim_v4, label='HII-CHI-mistry v4.20 grid', alpha=0.5)
    ax.scatter(x_HIIchim_v5, y_HIIchim_v5, label='HII-CHI-mistry v5.22 grid', alpha=0.5)

    ax.legend()
    ax.update({'xlabel': f'{line_labels[param]}', 'ylabel': line_labels[line]})
    plt.show()

param_list = ['logOH']
line_list = ['O3_5007A']

# # Plot the intervals
# for param, line in zip(param_list, line_list):
#     fig, ax = plt.subplots()
#     x_csv, y_csv = slice_grid(param, line, large_DF, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='O')
#     x_csv_N, y_csv_N = slice_grid(param, line, large_DF, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='N')
#     x_HIIchim_v5, y_HIIchim_v5 = slice_grid(param, line, HIICHIm_DF_v5, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
#     x_HIIchim_v2, y_HIIchim_v2 = slice_grid(param, line, HIICHIm_DF_v2, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
#
#     # ax.scatter(x_csv, (1-y_csv/y_HIIchim[1:])*100, label='CSV grid Carbon = O')
#     # ax.scatter(x_csv_N, (1-y_csv_N/y_HIIchim[1:])*100, label='CSV grid Carbon = N')
#
#     ax.legend()
#     ax.update({'xlabel': f'{line_labels[param]}', 'ylabel': f'{line_labels[line]} flux discrepancy % \n  with HII-CHI-mistry'})
#     plt.show()

# # [NII] discrepancies
# param_list = ['logNO', 'logNO']
# line_list = ['N2_6584A', 'O2_3726A_b']
# for param, line in zip(param_list, line_list):
#     fig, ax = plt.subplots()
#     x_csv, y_csv = slice_grid(param, line, large_DF, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='O')
#     x_csv_N, y_csv_N = slice_grid(param, line, large_DF, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='N')
#     x_HIIchim, y_HIIchim = slice_grid(param, line, HIICHIm_DF_v5, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)
#
#     ax.scatter(x_csv, (y_csv/y_HIIchim-1)*100, label='CSV grid Carbon = O')
#     ax.scatter(x_csv_N, (y_csv_N/y_HIIchim-1)*100, label='CSV grid Carbon = N')
#
#     ax.legend()
#     ax.update({'xlabel': f'{line_labels[param]}', 'ylabel': f'{line_labels[line]} flux discrepancy % \n  with HII-CHI-mistry'})
#     plt.show()

# NII fluxes versus N/O for various O
param_list = ['logNO'z]
line_list = ['N2_6584A', 'N2_6584A', 'O2_3726A_b', 'O2_3726A_b']

for param, line in zip(param_list, line_list):

    fig, ax = plt.subplots()

    for OH in [7.6, 7.8, 8.0]:

        x_csv, y_csv = slice_grid(param, line, large_DF, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='O')
        x_csv_N, y_csv_N = slice_grid(param, line, large_DF, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref, carbon='N')
        x_HIIchim_v5, y_HIIchim_v5 = slice_grid(param, line, HIICHIm_DF_v5, logOH=logOH_ref, logU=logU_ref, logNO=logNO_ref)


        ax.scatter(x_csv, y_csv, label='CSV grid Carbon = O', marker='*')
        ax.scatter(x_csv_N, y_csv_N, label='CSV grid Carbon = N', marker='*')
        ax.scatter(x_HIIchim_v2, y_HIIchim_v2, label='HII-CHI-mistry v2.00 grid', alpha=0.5)
        ax.scatter(x_HIIchim_v3, y_HIIchim_v3, label='HII-CHI-mistry v3.00 grid', alpha=0.5)
        ax.scatter(x_HIIchim_v4, y_HIIchim_v4, label='HII-CHI-mistry v4.20 grid', alpha=0.5)
        ax.scatter(x_HIIchim_v5, y_HIIchim_v5, label='HII-CHI-mistry v5.22 grid', alpha=0.5)

        ax.legend()
        ax.update({'xlabel': f'{line_labels[param]}', 'ylabel': line_labels[line]})
        plt.show()

param_list = ['logOH']
line_list = ['O3_5007A']