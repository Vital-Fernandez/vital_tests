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

# Configuration file
cfg_file = '../obsConf.ini'
obs_conf = lm.load_cfg(Path(cfg_file))
reduc_cfg = obs_conf['Megara_reduction']

# Data location
reduction_folder = Path(reduc_cfg['root_folder'])
data_folder = reduction_folder/'data'
rd_df_address = Path(reduc_cfg['rd_df_address'])

# Dataframe with files list
files_DF = lm.load_lines_log(f'{rd_df_address}')

# Targets lists
obj_list = reduc_cfg['obj_list']
std_list = reduc_cfg['std_star_list']

# Generate the task files for each OB:
OB_list = files_DF['OB'].unique()
OB_list.sort()
idx_start = 5

# Run the pipeline one OB and VPH at a time:
for OB in OB_list:

    if OB == 'OB0005':

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


                with open(req_yml) as f:
                    dataMap = yaml.safe_load(f)

                # Adjust manually the offsets
                for i, task_name in enumerate(task_list):
                    dataMap['database']['oblocks'][f'{task_list[i]}']['requirements']['extraction_offset'] = [1.5]
                with open(req_yml, 'w') as f:
                    yaml.dump(dataMap, f, sort_keys=False)

                # Define data manager
                dm = create_datamanager(req_yml, reduction_folder, data_folder)



                # req_yml = yaml.safe_load(str(req_yml))

                # Load the observation files
                with ctx.working_directory(reduction_folder):
                    sessions, loaded_obs = load_observations(task_file_list, is_session=False)
                    dm.backend.add_obs(loaded_obs)

                # Run the tasks
                for run_id in task_list:
                    print(f'Running: {run_id}')
                    output_run = run_reduce(dm, run_id)

