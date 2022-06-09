import numpy as np
import lime as lm
import pandas as pd
from megaradrp.tools.overplot_traces import main as megaradrp_traces_plot
from matplotlib import pyplot as plt
from pathlib import Path


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

        # Loop through the files types witht the extraction offset task parameter
        if master_traces_file.is_file():

            print(f'\n- File found for {OB}, {VPH}', master_traces_file)

            for file_type in ['arc']: #['arc', 'flat', 'object', 'stds']:

                # Find the files
                idcs_files = indexing_the_frames(files_DF, OB, VPH, file_type, bad_file_list)
                list_addresses = files_DF.loc[idcs_files].old_address.values # TODO old address in table is actually the new one

                # Start plotting the traces over the instrument frames if available
                if len(list_addresses) > 0:

                    for i, fits_address in enumerate(list_addresses):

                        # Confirm the file is there
                        fits_address = Path(fits_address)
                        if fits_address.is_file():

                            # Run the megardrp function
                            print(f'Treating {OB}-{VPH}: {file_type} ({fits_address.name})')
                            function_argument_list = [fits_address.as_posix(), master_traces_file.as_posix(), '--fibids']
                            function_argument_list += ['--bbox', '1900,2300,150,200']
                            megaradrp_traces_plot(function_argument_list)


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