import numpy as np
import pyneb as pn
import pandas as pd
import exoplanet as xo
import scipy as spy


def read_multidimensional_txt(file_address, file_headers, O3_label):
    # Read file
    grid_array = np.loadtxt(file_address, skiprows=1)
    lineLabels = file_headers[3:]

    # Recover axes dimensions
    axes_dimen = {}
    for i, ax_label in enumerate(file_headers[:3]):
        axes_dimen[ax_label] = np.unique(grid_array[:, i])

    # Sort the array
    idcs_sorted_grid = np.lexsort((grid_array[:, 2], grid_array[:, 1], grid_array[:, 0]))
    sorted_grid = grid_array[idcs_sorted_grid]

    # Rejoin per temp and density
    # for i, abund in enumerate(axes_dimen[O3_label]):
    grid_dict = {}
    for i, lineLabel in enumerate(lineLabels):
        grid_dict[lineLabel] = np.zeros((axes_dimen[file_headers[1]].size,
                                         axes_dimen[file_headers[2]].size,
                                         axes_dimen[file_headers[0]].size))

        for j, abund in enumerate(axes_dimen[O3_label]):

            idcsSubGrid = sorted_grid[:, 0] == abund
            lineGrid = sorted_grid[idcsSubGrid, 3:]
            # print(sorted_grid[idcsSubGrid, :3])
            lineCube = lineGrid.reshape((10, 9, lineLabels.size))
            grid_dict[lineLabel][:, :, j] = lineCube[:, :, i]

    # # Pre_analysis the recovering data:
    # idx_Line = 0
    # lineLabel = lineLabels[idx_Line]
    # cord_Te0_ne0 = grid_array[90, idx_Line + 3]
    # cord_Te1_ne0 = grid_array[91, idx_Line + 3]
    # cord_Te0_ne1 = grid_array[100, idx_Line + 3]
    # cord_Te1_ne1 = grid_array[101, idx_Line + 3]
    #
    # idx_abund = 1
    # print(f'cord 0 - 0 ({cord_Te0_ne0}): {grid_dict[lineLabel][0, 0, idx_abund] == cord_Te0_ne0} ({grid_dict[lineLabel][0, 0, idx_abund]})')
    # print(f'cord 1 - 0 ({cord_Te1_ne0}): {grid_dict[lineLabel][0, 1, idx_abund] == cord_Te0_ne1} ({grid_dict[lineLabel][0, 1, idx_abund]})')
    # print(f'cord 0 - 1 ({cord_Te0_ne1}): {grid_dict[lineLabel][1, 0, idx_abund] == cord_Te1_ne0} ({grid_dict[lineLabel][1, 0, idx_abund]})')
    # print(f'cord 1 - 1 ({cord_Te1_ne1}): {grid_dict[lineLabel][1, 1, idx_abund] == cord_Te1_ne1} ({grid_dict[lineLabel][1, 1, idx_abund]})')

    return grid_dict


# ---------------------------- Creating text grid file -----------------------------

# O3 = pn.Atom('O', 3)
# norm = 1e-20
#
# lineList = ['O3_4363A', 'O3_4959A', 'O3_5007A']
# lineWaves = np.array([4363.0, 4959.0, 5007.0])
#
# abund_array = np.arange(7.0, 8.0, 0.1)
# temp_array = np.arange(10000.0, 20000.0, 1000.0)
# den_array = np.arange(100, 1000, 100)
#
# grid_df = pd.DataFrame(columns=['Abund', 'Te', 'ne'] + lineList)
# grid_cube_orig = np.zeros((temp_array.size, den_array.size, abund_array.size))
#
# counter = 0
# for m, abund in enumerate(abund_array):
#     for n, den in enumerate(den_array):
#         for k, temp in enumerate(temp_array):
#
#             grid_df.loc[counter, 'Abund'] = abund
#             grid_df.loc[counter, 'Te'] = temp
#             grid_df.loc[counter, 'ne'] = den
#
#             for i, label in enumerate(lineList):
#                 emisLine = abund * (O3.getEmissivity(tem=temp, den=den, wave=lineWaves[i])/norm)
#                 grid_df.loc[counter, lineList[i]] = emisLine
#                 if i == 0:
#                     grid_cube_orig[k, n, m] = emisLine
#
#             counter += 1
#
# matrix_array = abund_array[0] * (O3.getEmissivity(tem=temp_array, den=den_array, wave=lineWaves[0])/norm)
# matrix_array_00 = abund_array[0] * (O3.getEmissivity(tem=temp_array[0], den=den_array[0], wave=lineWaves[0])/norm)
# matrix_array_01 = abund_array[0] * (O3.getEmissivity(tem=temp_array[0], den=den_array[1], wave=lineWaves[0])/norm)
# output_file = 'D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/treatment/grid_file.txt'

# Safe the dataframe
# with open(output_file, 'wb') as txt_file:
#     string_DF = grid_df.to_string(index=False)
#     txt_file.write(string_DF.encode('UTF-8'))


# ---------------------------- Recovering grid file -----------------------------

O3 = pn.Atom('O', 3)
norm = 1e-20
true_values = np.array([18000.0, 400.0, 7.5])

lineList = ['O3_4363A', 'O3_4959A', 'O3_5007A']
lineWaves = np.array([4363.0, 4959.0, 5007.0])

abund_array = np.array([7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9])
temp_array = np.arange(10000.0, 20000.0, 1000.0)
den_array = np.arange(100, 1000, 100)

# Compute the cube
emisValues = O3.getEmissivity(temp_array, den_array, wave=4363.0, product=True)
emisCube = np.ones((temp_array.size, den_array.size, abund_array.size))
for i, abund_value in enumerate(abund_array):
    emisCube[:, :, i] = (emisValues/norm) * abund_value

# Get the data
output_file = 'D:/Dropbox/Astrophysics/Papers/gtc_greenpeas/treatment/grid_file.txt'
headers = np.array(['Abund', 'Te', 'ne', 'O3_4363A', 'O3_4959A', 'O3_5007A'])
grid_dict = read_multidimensional_txt(output_file, headers, O3_label='Abund')

# Compare calculations
emis_true = true_values[2] * O3.getEmissivity(tem=true_values[0], den=true_values[1], wave=4363.0) / norm
print('1 True value', emis_true)

i_temp, i_den, i_abund = temp_array == true_values[0], den_array == true_values[1], abund_array == true_values[2]
emis_grid_1 = grid_dict['O3_4363A'][i_temp, i_den, i_abund]
emis_grid_2 = emisCube[i_temp, i_den, i_abund]

print('2 Manual grid file', emis_grid_1)
print('3 Manual grid array', emis_grid_2)

# Scipy interpolation RegularGridInterpolator
spy_RGridInterp_file = spy.interpolate.RegularGridInterpolator((temp_array, den_array, abund_array), grid_dict['O3_4363A'])
spy_RGridInterp_array = spy.interpolate.RegularGridInterpolator((temp_array, den_array, abund_array), emisCube)
print('4 Scipy RegularGridInterpolator file', spy_RGridInterp_file([[18000.0, 400.0, 7.5]]))
print('5 Scipy RegularGridInterpolator array', spy_RGridInterp_array([[18000.0, 400.0, 7.5]]))

exop_interp = xo.interp.RegularGridInterpolator([temp_array, den_array, abund_array], emisCube, nout=1)
coordB = np.stack(([18000.0], [400.0], [7.5]), axis=-1)
print('6 Exoplanet interpolation array', exop_interp.evaluate(coordB).eval())


# idcsCubeGridt = np.lexsort((grid_array[:, 0], grid_array[:, 1], grid_array[:, 2]))
# gridSorted = grid_array[idcsCubeGridt]
# for row in gridSorted:
#     print(row[0], row[1], row[2])
#

# Pre_analysis the recovering data:
# idx_Line = 0
# cord_Te0_ne0 = grid_array[0, idx_Line + 3]
# cord_Te1_ne0 = grid_array[1, idx_Line + 3]
# cord_Te0_ne1 = grid_array[10, idx_Line + 3]
# cord_Te1_ne1 = grid_array[11, idx_Line + 3]
#
# # Reordering the text file to high dimensional array
# print(f'\nabund = {abund_array[0]}')
# idcs_sorted_grid = np.lexsort((grid_array[:, 2], grid_array[:, 1], grid_array[:, 0]))
# sorted_grid = grid_array[idcs_sorted_grid]
# idcsSubGrid = sorted_grid[:, 0] == abund_array[0]
# lineGrid = sorted_grid[idcsSubGrid, 3:]
# print(sorted_grid[idcsSubGrid, :3])
# lineCube = lineGrid.reshape((10, 9, 3))
#
# print(f'cord 0 - 0 ({cord_Te0_ne0}): {lineCube[0, 0, idx_Line] == cord_Te0_ne0} ({lineCube[0, 0, idx_Line]})')
# print(f'cord 1 - 0 ({cord_Te1_ne0}): {lineCube[0, 1, idx_Line] == cord_Te0_ne1} ({lineCube[0, 1, idx_Line]})')
# print(f'cord 0 - 1 ({cord_Te0_ne1}): {lineCube[1, 0, idx_Line] == cord_Te1_ne0} ({lineCube[1, 0, idx_Line]})')
# print(f'cord 1 - 1 ({cord_Te1_ne1}): {lineCube[1, 1, idx_Line] == cord_Te1_ne1} ({lineCube[1, 1, idx_Line]})')

# Reordering the text file to high dimensional array
# idcs_sorted_grid = np.lexsort((grid_array[:, 2], grid_array[:, 0], grid_array[:, 1]))
# sorted_grid = grid_array[idcs_sorted_grid]
# lineGrid = sorted_grid[:, 3:]
#
# line_4Darray = lineGrid.reshape((9, 10, 10, 3))
# lineCubeAbund0 = line_4Darray[:, :, 0, idx_Line]
#
# print(f'\nTemperature line {lineCubeAbund0[0, :]}')
# print(f'First temps {cord_Te0_ne0, cord_Te1_ne0}\n')
# print(f'Density line {lineCubeAbund0[:, 0]}')
# print(f'first dens {cord_Te0_ne0, cord_Te0_ne1}\n')
#
# print(f'cord 0 - 0 ({cord_Te0_ne0}): {lineCubeAbund0[0, 0] == cord_Te0_ne0} ({lineCubeAbund0[0, 0]})')
# print(f'cord 1 - 0 ({cord_Te1_ne0}): {lineCubeAbund0[0, 1] == cord_Te0_ne1} ({lineCubeAbund0[0, 1]})')
# print(f'cord 0 - 1 ({cord_Te0_ne1}): {lineCubeAbund0[1, 0] == cord_Te1_ne0} ({lineCubeAbund0[1, 0]})')
# print(f'cord 1 - 1 ({cord_Te1_ne1}): {lineCubeAbund0[1, 1] == cord_Te1_ne1} ({lineCubeAbund0[1, 1]})')


# for OH in grid_axes['abund']:
#     print(f'\nabund = {OH}')
#     idcsSubGrid = grid_array[:, 0] == OH
#     subGrid = grid_array[idcsSubGrid, :]
#     lineGrid = grid_array[idcsSubGrid, 3:]
#     for Temp in np.unique(subGrid[:, 1]):
#         idcs_Teff = subGrid[:, 1] == Temp
#         ne_array = np.unique(subGrid[idcs_Teff, 2])
#         print(f'Given Teff = {Temp} -> ({ne_array.size}) logU array ({np.min(ne_array)}, {np.max(ne_array)})')
