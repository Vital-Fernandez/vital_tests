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


def delete_task_temp_folder(ob_id, task_id, root_dir):

    for folder in ['work', 'result']:
        temp_folder = root_dir/f'obsid{ob_id}_{task_id}_{folder}'
        if temp_folder.is_dir():
            shutil.rmtree(temp_folder)

    return

# Loading project configuration file
obs_conf = lm.load_cfg(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\obsConf.ini')
reduction_cfg = obs_conf['Megara_reduction']
obj_list = reduction_cfg['obj_list']
std_list = reduction_cfg['std_star_list']

# Dataframe with files list
rd_df_address = Path(reduction_cfg['rd_df_address'])
files_DF = lm.load_lines_log(f'{rd_df_address}.txt')

# Stating folder structure
instructions_folder = Path(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\MEGARA_reduction')
reduction_folder = Path(reduction_cfg['root_folder'])
data_folder = reduction_folder/'data'

# Generate the task files for each OB:
OB_list = files_DF['OB'].unique()
OB_list.sort()
idx_start = 6

# Run the pipeline one OB and VPH at a time:
for OB in OB_list:

    idcs_OB = files_DF.OB == OB
    VPH_list = files_DF.loc[idcs_OB, 'VPH'].unique()

    for VPH in VPH_list:
        if VPH != 'MR-B':

            task_file_address = f'{reduction_folder}/{OB}_{VPH}_task_list.txt'
            task_DF = pd.read_csv(task_file_address, delim_whitespace=True, header=0, index_col=0)

            # Create clean requirements file at each run
            input_yml = reduction_folder/f'control_{OB}_{VPH}_phase1.yaml'
            req_yml = reduction_folder/f'control_{OB}_{VPH}_phase2.yaml'

            # Copy previous yml so the original is not rewritten in current phase
            shutil.copyfile(input_yml, req_yml)

            # Decide tasks to run
            idcs_tasks = task_DF.index >= idx_start
            task_list = task_DF.loc[idcs_tasks].task_id.values
            task_file_list = task_DF.loc[idcs_tasks].file_name.values
            idcs_tasks = task_DF.loc[idcs_tasks].index.values

            # Define data manager
            dm = create_datamanager(req_yml, reduction_folder, data_folder)

            with open(req_yml) as f:
                # use safe_load instead load
                dataMap = yaml.safe_load(f)

            req_yml = yaml.safe_load(str(req_yml))

            # Load the observation files
            with ctx.working_directory(reduction_folder):
                sessions, loaded_obs = load_observations(task_file_list, is_session=False)
                dm.backend.add_obs(loaded_obs)

            # Run the tasks
            for run_id in task_list:
                print(f'Running: {run_id}')
                output_run = run_reduce(dm, run_id)

