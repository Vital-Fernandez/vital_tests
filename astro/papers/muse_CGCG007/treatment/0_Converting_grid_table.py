import numpy as np
import pandas as pd
import lime
from pathlib import Path
import src.specsiser as sr

columns_renam_dict = {'12+log(O/H)': 'logOH',
                      'log(N/O)': 'logNO',
                      'log(U)': 'logU',
                      'OII3726A': 'O2_3726A',
                      'OII3729A': 'O2_3729A',
                      'NeIII3869A': 'Ne3_3869A',
                      'OIII4363A': 'O3_4363A',
                      'HeI4471A': 'He1_4471A',
                      'HeII4686A': 'He2_4686A',
                      'ArIV4740A': 'Ar4_4740A',
                      'OIII4959A': 'O3_4959A',
                      'OIII5007A': 'O3_5007A',
                      'NII5755A': 'N1_5755A',
                      'HeI5876A': 'HeI_5876A',
                      'SIII6312A': 'S3_6312A',
                      'NII6548A': 'N2_6548A',
                      'HI6563A': 'H1_6563A',
                      'NII6583A': 'N2_6584A',
                      'HeI6678A': 'He1_6678A',
                      'SII6716A': 'S2_6716A',
                      'SII6731A': 'S2_6731A',
                      'ArII7135A': 'Ar3_7136A',
                      'OII7319A': 'O2_7319A',
                      'OII7330A': 'O2_7330A',
                      'ArII7751A': 'Ar3_7751A',
                      'SIII9069A': 'S3_9069A',
                      'SIII9532A': 'S3_9531A'}

orig_grid_file = Path(r'D:\Dropbox\Astrophysics\Papers\muse_CGCG007\data\C17_Popstar_1Myr.dat')
output_grid_file = Path(r'D:\Dropbox\Astrophysics\Papers\muse_CGCG007\data\formated_log_C17_Popstar_1Myr.dat')

grid = pd.read_csv(orig_grid_file, delim_whitespace=True)
grid.rename(columns=columns_renam_dict, inplace=True)
columns = list(grid.columns)

# Replacing the log(N/O) = 1.7750 for 1.8750 at 12+log(O/H) = 7.5
# Some logOH = 7.6000000 logN/O = -2.0000000 entries are repeated

idcs_wrong_NO = (grid.logOH == 7.5) & (grid.logNO == -1.7750)
grid.loc[idcs_wrong_NO, 'logNO'] = -1.8750


# Remove Unecesary entries
logOH_good = [7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.,  8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9., 9.1]
logNO_good = [-2., -1.875, -1.75, -1.625, -1.5, -1.375, -1.25, -1.125, -1., -0.875, -0.75,  -0.625, -0.5, -0.375, -0.25, -0.125, 0.]

idcs_good = grid.logOH.isin(logOH_good) & grid.logNO.isin(logNO_good)
grid = grid.loc[idcs_good]

for OH in logOH_good:
    idcs_good = grid.logOH == OH
    slice_grid = grid.loc[idcs_good]
    # print(f'{OH} -> {np.sum(idcs_good)} entries')
    print(f'{OH} -> {np.unique(slice_grid.logNO).size} log(NO); {np.unique(slice_grid.logU).size} log(U)')

# Move the third column to the first
grid = grid[['logOH', 'logU', 'logNO'] + columns[3:]]

# Sort by "logOH", "logU", "logNO" axes
grid.sort_values(["logOH", "logU", "logNO"], ascending=True, inplace=True)

# Convert to loc scale
grid.iloc[:, 3:] = np.log10(grid.iloc[:, 3:])
grid.replace(-np.inf, -6.0000, inplace=True)

print(grid)
print(len(grid))

gw = sr.GridWrapper()
grid_dict, axes_cords_a = gw.ndarray_from_DF(grid, axes_columns=["logOH", "logU", "logNO"])

# Safe the dataframe
with open(output_grid_file, 'wb') as txt_file:
    string_DF = grid.to_string(index=False)
    txt_file.write(string_DF.encode('UTF-8'))
