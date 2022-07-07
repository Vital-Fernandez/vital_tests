import numpy as np
import pandas as pd
import lime as lm
import shutil
import yaml

from pathlib import Path
from numina.user.helpers import create_datamanager, load_observations
from numina.user.baserun import run_reduce
from numina.util import context as ctx
from astropy.io import fits
from timeit import default_timer as timer
from .tools import indexing_the_frames, store_task

# Configuration file
cfg_file = '../obsConf.ini'
obs_conf = lm.load_cfg(Path(cfg_file))
reduc_cfg = obs_conf['Megara_reduction']
task_list = reduc_cfg['task_list']

# Data location
reduction_folder = Path(reduc_cfg['root_folder'])
rd_df_address = Path(reduc_cfg['rd_df_address'])
data_folder = reduction_folder/'data'
original_yml = Path().resolve()/'shoc579_req.yml'

# Dataframe with files list
files_DF = lm.load_lines_log(f'{rd_df_address}')
bad_file_list = np.loadtxt(reduc_cfg['issue_frames_file'], usecols=0, skiprows=1, dtype=str)

# Generate the task files for each OB:
OB_list = files_DF['OB'].unique()
OB_list.sort()

# Parameters for the current task
task = 'bias'
target_list = [None]
i_task = task_list.index(task)

# Run the pipeline one OB and VPH at a time:
start = timer()
for OB in OB_list:

    # The VPHs in this OB but exclude the MR-B which are only used in the Bias operation
    idcs_OB = files_DF.OB == OB
    VPH_list = files_DF.loc[idcs_OB, 'VPH'].unique()
    idx_MR_B = np.where(VPH_list == 'MR-B')[0][0]
    VPH_list = np.delete(VPH_list, idx_MR_B)

    # Loop through the VPHs
    for VPH in VPH_list:

        # Loop through the sky targets (if applicable)
        for target in target_list:

            # Task label and number
            task_ID = f'{OB}_{VPH}_{task}' if target is None else f'{OB}_{VPH}_{task}_{target}'

            # Get fits files for the operation
            VPH_in = 'MR-B' if task == 'bias' else VPH
            task_in = 'flat' if task in ['trace_map', 'model_map', 'fiber_flat'] else task
            idcs_fits = indexing_the_frames(files_DF, OB, VPH_in, task_in, target, bad_file_list)

            # Requirements file making a copy from the previous step
            req_yml = f'phase{i_task}_{task_ID}'

            # Task configuration
            task_conf = None

            # Store the task configuration into a configuration file if there are files for it in the OB
            if np.sum(idcs_fits):
                req_dict = {} if task_conf is None else {'requirements': task_conf}
                store_task(task_ID, VPH, reduction_folder, files_DF, task, req_dict, idcs_fits, reduc_cfg[f'{VPH}_lamp'])

        # Create copy of previous task requirements ymal for the current operation
        current_yml = reduction_folder/f'phase{i_task}_{task}_{OB}_{VPH}.yaml'
        previous_yml = current_yml if i_task == 0 else reduction_folder/f'phase{i_task}_{task}_{OB}_{VPH}.yaml'

        # If the task was run before delete previous results


        # Define data manager
        dm = create_datamanager(reduction_folder/current_yml, reduction_folder, data_folder)

        # Load the observation files
        with ctx.working_directory(reduction_folder):
            sessions, loaded_obs = load_observations(task_file_list, is_session=False)
            dm.backend.add_obs(loaded_obs)

        # # Run the tasks
        # for i, idx_task in enumerate(idcs_tasks):
        #     if idx_task <= idx_finish:
        #         run_id = task_list[i]
        #         print(f'\n=================================Running: {run_id}================================\n')
        #         warning_messange(run_id, reduction_folder)
        #         run_reduce(dm, run_id)
        #
end = timer()

print(f'Working time: {(end-start)/60:0.1f} mins')
