import pandas as pd
import lime as lm
import yaml
from astro.papers.SHOC579_project.SHOC579_methods import megaradrp_modes
from pathlib import Path


def store_task(task_ID_ref, VPH_ref, output_folder, DF_files, task_name, DF_tasks, idx_task, conf_task):

    # Generate task yml
    yml_dict = {'id': task_ID_ref,
                'mode': megaradrp_modes[task_name],
                'instrument': 'MEGARA',
                'frames': list(DF_files.loc[idcs_fits].index.values)}
    yml_dict.update(conf_task)

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
reduction_cfg = obs_conf['Megara_reduction']

# Data location
reduction_folder = Path(reduction_cfg['root_folder'])
data_folder = reduction_folder/'data'
rd_df_address = Path(reduction_cfg['rd_df_address'])

# Loading project configuration file
obj_list = reduction_cfg['obj_list']
std_list = reduction_cfg['std_star_list']
bad_file_list = reduction_cfg['issue_frames_list']

# Dataframe with files list
files_DF = lm.load_lines_log(f'{rd_df_address}.txt')

# Generate the task files for each OB:
OB_list = files_DF['OB'].unique()
OB_list.sort()

# List of tasks in the reduction
# task_list = ['bias', 'trace_map', 'arc', 'fiber_flat', 'lcb_std', 'lcb_image']
task_list = ['bias', 'trace_map', 'model_map', 'arc', 'fiber_flat', 'lcb_std', 'lcb_image']

for OB in OB_list:

    idcs_OB = files_DF.OB == OB
    VPH_list = files_DF.loc[idcs_OB, 'VPH'].unique()

    for VPH in VPH_list:

        # Exclude the MR-B
        if VPH != 'MR-B':

            # Dataframe for storing the tasks:
            i_task = 1
            task_DF = pd.DataFrame(columns=['task_id', 'file_name', 'VPH'])

            for task in task_list:
                extra_conf = {}

                # File selection
                if task == 'bias':
                    task_ID = f'{OB}_{VPH}_{i_task}_{task}'
                    idcs_fits = (files_DF.OB == OB) & (files_DF.VPH == 'MR-B') & (files_DF.type == task) \
                                & (files_DF.reduc_tag == 'raw') & (~files_DF.index.isin(bad_file_list))

                    store_task(task_ID, VPH, reduction_folder, files_DF, task, task_DF, i_task, extra_conf)
                    i_task += 1

                if task == 'arc':
                    task_ID = f'{OB}_{VPH}_{i_task}_{task}'
                    idcs_fits = (files_DF.OB == OB) & (files_DF.VPH == VPH) & (files_DF.type == task) \
                                & (files_DF.reduc_tag == 'raw') & (~files_DF.index.isin(bad_file_list))

                    extra_conf['requirements'] = {'extraction_offset':  [0.0],
                                                  'store_pdf_with_refined_fits': 1}

                    store_task(task_ID, VPH, reduction_folder, files_DF, task, task_DF, i_task, extra_conf)
                    i_task += 1

                if task in ['trace_map', 'model_map', 'fiber_flat']:
                    task_ID = f'{OB}_{VPH}_{i_task}_{task}'
                    idcs_fits = (files_DF.OB == OB) & (files_DF.VPH == VPH) & (files_DF.type == 'flat') \
                                & (files_DF.reduc_tag == 'raw') & (~files_DF.index.isin(bad_file_list))

                    if task == 'fiber_flat':
                        extra_conf['requirements'] = {'extraction_offset': [0.0]}

                    store_task(task_ID, VPH, reduction_folder, files_DF, task, task_DF, i_task, extra_conf)
                    i_task += 1

                if task == 'lcb_std':
                    for std in std_list:
                        if files_DF.loc[idcs_OB, 'object'].str.contains(std).any():

                            task_ID = f'{OB}_{VPH}_{i_task}_{std}_{task}'

                            idcs_fits = (files_DF.OB == OB) & (files_DF.VPH == VPH) & \
                                        (files_DF.object == std) & (files_DF.reduc_tag == 'raw') \
                                        & (~files_DF.index.isin(bad_file_list))

                            req_dict = {'extraction_offset': [0.0],
                                        'reference_extinction': 'extinction_LP.txt'}

                            if std == 'HR8634':
                                req_dict['reference_spectrum'] = 'mhr8634.dat'

                            if std == 'HR7596':
                                req_dict['reference_spectrum'] = 'mhr7596.dat'

                            extra_conf['requirements'] = req_dict

                            store_task(task_ID, VPH, reduction_folder, files_DF, task, task_DF, i_task, extra_conf)
                            i_task += 1

                if task == 'lcb_image':
                    for obj in obj_list:
                        if files_DF.loc[idcs_OB, 'object'].str.contains(obj).any():

                            task_ID = f'{OB}_{VPH}_{i_task}_{obj}_{task}'

                            idcs_fits = (files_DF.OB == OB) & (files_DF.VPH == VPH) & \
                                        (files_DF.object == obj) & (files_DF.reduc_tag == 'raw') \
                                        & (~files_DF.index.isin(bad_file_list))

                            extra_conf['requirements'] = {'extraction_offset': [0.0],
                                                          'reference_extinction': 'extinction_LP.txt'}

                            store_task(task_ID, VPH, reduction_folder, files_DF, task, task_DF, i_task, extra_conf)
                            i_task += 1


            # Save task DF to a text file
            with open(f'{reduction_folder}\{OB}_{VPH}_task_list.txt', 'wb') as output_file:
                string_DF = task_DF.to_string()
                output_file.write(string_DF.encode('UTF-8'))