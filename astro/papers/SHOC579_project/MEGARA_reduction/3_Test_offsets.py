import numpy as np
import lime as lm
import pandas as pd
import shutil
import os

from megaradrp.tools.overplot_traces import main as megaradrp_traces_plot
from matplotlib import pyplot as plt, rcParams
from pathlib import Path


rcParams.update({'figure.figsize': (14, 7)})

def indexing_the_frames(log, OB_task, VPH_task, type_task, exclude_list=[]):

    idcs_frames = (log.OB == OB_task) & (log.VPH == VPH_task) & (log.type == type_task) & (~log.index.isin(exclude_list))

    return idcs_frames

# Configuration file
cfg_file = '../obsConf.ini'
obs_conf = lm.load_cfg(Path(cfg_file))
reduction_cfg = obs_conf['Megara_reduction']

# Data location
reduction_folder = Path(reduction_cfg['root_folder'])
data_folder = reduction_folder/'data'
rd_df_address = Path(reduction_cfg['rd_df_address'])
bad_file_list = np.loadtxt(reduction_cfg['issue_frames_file'], usecols=0, skiprows=1, dtype=str)

# Dataframe with files list
files_DF = pd.read_csv(rd_df_address, delim_whitespace=True, header=0, index_col=0) #lm.load_lines_log(f'{rd_df_address}')
OB_list = files_DF['OB'].unique()
OB_list.sort()

# Run the pipeline one OB and VPH at a time:
counter = 0
for OB in OB_list:

    #Get list of VPH in OB
    idcs_OB = files_DF.OB == OB
    VPH_list = files_DF.loc[idcs_OB, 'VPH'].unique()

    # Exclude the MR-B which are only used in the Bias operation
    idx_MR_B = np.where(VPH_list == 'MR-B')[0][0]
    VPH_list = np.delete(VPH_list, idx_MR_B)

    # Loop through the VPH to get the master traces file (if available)
    for VPH in VPH_list:

        task_file_address = f'{reduction_folder}/{OB}_{VPH}_task_list.txt'
        task_DF = pd.read_csv(task_file_address, delim_whitespace=True, header=0, index_col=0)
        master_traces_file = Path(f'{reduction_folder}/obsid{OB}_{VPH}_trace_map_2_result/master_traces.json')
        master_traces_ds9_file = Path(f'{reduction_folder}/obsid{OB}_{VPH}_model_map_3_work/ds9_raw.reg')

        # Loop through the files types witht the extraction offset task parameter
        if master_traces_file.is_file():

            # print(f'\n- File found for {OB}, {VPH}', master_traces_file)

            for file_type in ['flat', 'arc', 'object', 'stds']:

                print(f'\n-------------------------------- ({counter}) Treating {OB}-{VPH}: {file_type} --------------------------\n')

                # Find the files
                idcs_files = indexing_the_frames(files_DF, OB, VPH, file_type, bad_file_list)
                list_addresses = files_DF.loc[idcs_files].address.values # TODO old address in table is actually the new one

                # Start plotting the traces over the instrument frames if available
                if len(list_addresses) > 0:

                    for i, fits_address in enumerate(list_addresses):

                        # Confirm the file is there
                        fits_address = Path(fits_address)
                        if fits_address.is_file():

                            # Check if there is an offset for this file
                            stored_offset = obs_conf['Megara_reduction'].get(f'{OB}_{VPH}_{file_type}', None)

                            # Run the megardrp function
                            # print(f'-------------------------------- ({counter}) Treating {OB}-{VPH}: {file_type} --------------------------')
                            # print(f'{counter} {fits_address} ({i+1}/{len(list_addresses)})')
                            # print(f'{OB}_{VPH}_{file_type}_offset={stored_offset}\n')

                            # # Copy previous yml so the original is not rewritten in current phase
                            # copy_fits = reduction_folder/'test_extensions'/f'{file_type}'/fits_address.name
                            # copy_master = reduction_folder/'test_extensions'/f'{file_type}'/f'{OB}_{VPH}.reg'
                            #

                            # # shutil.copyfile(fits_address, copy_fits)
                            # # shutil.copyfile(master_traces_file, copy_master)
                            # os.system(ds9_command)

                            # '/mnt/AstroData/Observations/SHOC579/MEGARA/test_extensions'
                            #
                            # shutil.copyfile(input_yml, req_yml)

                            function_argument_list = [fits_address.as_posix(), master_traces_file.as_posix(), '--fibids', '--rawimage']
                            # function_argument_list += ['--bbox', '1850,2350,125,225']
                            megaradrp_traces_plot(function_argument_list)

                            ds9_command = f'ds9 -tile {fits_address} -cmap Heat -zscale -regions {master_traces_ds9_file}'
                            print(ds9_command)
                            os.system(ds9_command)

                            #
                            # function_argument_list += ['--bbox', '1900,2200,3850,4000']
                            # megaradrp_traces_plot(function_argument_list)
                            #
                            # counter += 1


# Reading apertures from trace file
#
# import numina.types.structured as structured
#
# apers = structured.open(trace_file_address)
# for geot in apers.contents:
#     if geot.valid:
#         fibid = geot.fibid                                      # Fiber ID
#         center_model = geot.aper_center()                       # Numpy Polynomial
#         num = int(float(4092 - 4 + 1) + 0.5)                    # Number of fiber plot elements
#         xp = np.linspace(start=4, stop=4092, num=num) + 51      # X values lines
#         yp = center_model(xp) + 1                               # y values lines
#         print(fibid)
#         if fibid == 1:
#             fig, ax = plt.subplots()
#             ax.plot(xp, yp)
#             ax.legend()
#             plt.show()