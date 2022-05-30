import numpy as np
import pandas as pd
import lime as lm
import yaml
from astro.papers.SHOC579_project.SHOC579_methods import megaradrp_modes
from pathlib import Path


def indexing_the_frames(log, OB_task, VPH_task, type_task, obj_task=None, exclude_list=[]):

    if obj_task is None:
        idcs_frames = (log.OB == OB_task) & (log.VPH == VPH_task) & (log.type == type_task) & (~log.index.isin(exclude_list))
    else:
        idcs_frames = (log.OB == OB_task) & (log.VPH == VPH_task) & (log.object == obj_task) & (~log.index.isin(exclude_list))

    return idcs_frames


def store_task(task_ID_ref, VPH_ref, output_folder, DF_files, task_name, DF_tasks, idx_task, conf_task, idcs_files):

    # Generate task yml
    yml_dict = {'id': task_ID_ref,
                'mode': megaradrp_modes[task_name],
                'instrument': 'MEGARA',
                'frames': list(DF_files.loc[idcs_files].index.values)}
    yml_dict.update(conf_task)

    print(f'\n --- {task_ID_ref} ({np.sum(idcs_files)} files)--- ')
    for file in DF_files.loc[idcs_files].index.values:
        print(file)

    # Save yml to a text file
    dict_adress = f'{output_folder}/{task_ID_ref}.yml'
    with open(dict_adress, 'w') as f:
        yaml.dump(yml_dict, f, sort_keys=False)

    # Store DataFrame data
    DF_tasks.loc[idx_task, 'task_id'] = task_ID_ref
    DF_tasks.loc[idx_task, 'file_name'] = f'{task_ID_ref}.yml'
    DF_tasks.loc[idx_task, 'VPH'] = VPH_ref

    return


# Configuration file
cfg_file = '../obsConf.ini'
obs_conf = lm.load_cfg(Path(cfg_file))
reduc_cfg = obs_conf['Megara_reduction']

# Data location
reduction_folder = Path(reduc_cfg['root_folder'])
data_folder = reduction_folder/'data'
rd_df_address = Path(reduc_cfg['rd_df_address'])

# Loading project configuration file
obj_list = reduc_cfg['obj_list']
std_list = reduc_cfg['std_star_list']
bad_file_list = np.loadtxt(reduc_cfg['issue_frames_file'], usecols=0, skiprows=1, dtype=str)

# Dataframe with files list
files_DF = lm.load_lines_log(f'{rd_df_address}')

# Generate the task files for each OB:
OB_list = files_DF['OB'].unique()
OB_list.sort()

# List of tasks in the reduction
task_list = ['bias', 'trace_map', 'model_map', 'arc', 'fiber_flat', 'lcb_acq', 'lcb_std', 'lcb_image']

for OB in OB_list:

    idcs_OB = files_DF.OB == OB
    VPH_list = files_DF.loc[idcs_OB, 'VPH'].unique()

    # Exclude the MR-B which are only used in the Bias operation
    idx_MR_B = np.where(VPH_list == 'MR-B')[0][0]
    VPH_list = np.delete(VPH_list, idx_MR_B)

    for VPH in VPH_list:

        # Dataframe for storing the tasks:
        i_task = 0
        task_DF = pd.DataFrame(columns=['task_id', 'file_name', 'VPH'])
        task_conf = None

        # Loop throught the tasks
        for task in task_list:

            if task in ['lcb_acq', 'lcb_std']:
                target_list = std_list
            elif task == 'lcb_image':
                target_list = obj_list
            else:
                target_list = [None]

            # Loop through the objects
            for target in target_list:

                # Task label
                task_ID = f'{OB}_{VPH}_{task}' if target is None else f'{OB}_{VPH}_{task}_{target}'

                # Get the files for each operation
                VPH_in = 'MR-B' if task == 'bias' else VPH
                task_in = 'flat' if task in ['trace_map', 'model_map', 'fiber_flat'] else task
                idcs_fits = indexing_the_frames(files_DF, OB, VPH_in, task_in, target, bad_file_list)

                # Get the task configuration

                # Arcs conf
                if task == 'arc':
                    task_conf = {'extraction_offset': reduc_cfg.get(f'{OB}_{VPH}_{task}_extraction_offset', [0.0]),
                                 'store_pdf_with_refined_fits': reduc_cfg.get(f'{OB}_{VPH}_{task}_extraction_offset', 1)}

                # Flats conf
                if task in ['trace_map', 'model_map', 'fiber_flat']:
                    task_conf = {'extraction_offset': reduc_cfg.get(f'{OB}_{VPH}_{task}_extraction_offset', [0.0])}

                # Flux calibration
                if task in ['lcb_acq', 'lcb_std', 'lcb_image']:
                    extraction_offset = reduc_cfg.get(f'{OB}_{VPH}_{task}_{target}_extraction_offset', [0.0])
                    ref_extincion = 'extinction_LP.txt'

                    # Configuration of the task
                    task_conf = {'extraction_offset': extraction_offset}

                    # Standard star
                    if task == 'lcb_std':
                        reference_spectrum = reduc_cfg[f'{target}_{task}_reference_spectrum']
                        task_conf['reference_extinction'] = ref_extincion
                        task_conf['reference_spectrum'] = reference_spectrum
                        task_conf['sigma_resolution'] = reduc_cfg.get(f'{OB}_{VPH}_{task}_{target}_extraction_offset', 50)

                    # Science target
                    if task == 'lcb_image':
                        task_conf['reference_extinction'] = ref_extincion

                # Store the task configuration into a configuration file if there are files for it in the OB
                if np.sum(idcs_fits):
                    req_dict = {} if task_conf is None else {'requirements': task_conf}
                    store_task(task_ID, VPH, reduction_folder, files_DF, task, task_DF, i_task, req_dict, idcs_fits)
                    i_task += 1

        # Save table with the task list per OB-VPH
        with open(f'{reduction_folder}/{OB}_{VPH}_task_list.txt', 'wb') as output_file:
            string_DF = task_DF.to_string()
            output_file.write(string_DF.encode('UTF-8'))
