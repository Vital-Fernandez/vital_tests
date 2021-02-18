import numpy as np

HCm_folder = 'D:/Dropbox/Astrophysics/Tools/HCm-Teff_v5.01/'
HCm_grid = 'C17_bb_Teff_30-90_pp.dat'

grid_array = np.loadtxt(HCm_folder+HCm_grid)

grid_axes = dict(OH=np.unique(grid_array[:, 0]),
              Teff=np.unique(grid_array[:, 1]),
              logU=np.unique(grid_array[:, 2]))

for ax_label, ax_values in grid_axes.items():
    print(f'Grid {ax_label}: {ax_values}')

# idcsCubeGridt = np.lexsort((grid_array[:, 0], grid_array[:, 1], grid_array[:, 2]))
# gridSorted = grid_array[idcsCubeGridt]
# for row in gridSorted:
#     print(row[0], row[1], row[2])
#
#
# for OH in grid_axes['OH']:
#     print(f'\nOH = {OH}')
#     idcsSubGrid = grid_array[:, 0] == OH
#     subGrid = grid_array[idcsSubGrid, :]
#     for Teff in np.unique(subGrid[:, 1]):
#         idcs_Teff = subGrid[:, 1] == Teff
#         logU_array = np.unique(subGrid[idcs_Teff, 2])
#         print(f'Given Teff = {Teff} -> ({logU_array.size}) logU array ({np.min(logU_array)}, {np.max(logU_array)})')
#
#




