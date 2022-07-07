import numpy as np
import shutil
import yaml
from astropy.io import fits


megaradrp_modes = {'bias'       : 'MegaraBiasImage',
                   'arc'        : 'MegaraArcCalibration',
                   'trace_map'  : 'MegaraTraceMap',
                   'slit_flat'  : 'MegaraSlitFlat',
                   'model_map'  : 'MegaraModelMap',
                   'fiber_flat' : 'MegaraFiberFlatImage',
                   'lcb_acq'    : 'MegaraLcbAcquisition',
                   'lcb_std'    : 'MegaraLcbStdStar',
                   'lcb_image'  : 'MegaraLcbImage'}

task_order = {'bias'       : 1,
              'arc'        : 2,
              'trace_map'  : 3,
              'slit_flat'  : 4,
              'model_map'  : 5,
              'fiber_flat' : 6,
              'lcb_acq'    : 7,
              'lcb_std'    : 8,
              'lcb_image'  : 9}


def delete_task_temp_folder(ob_id, task_id, root_dir):

    for folder in ['work', 'result']:
        temp_folder = root_dir/f'obsid{ob_id}_{task_id}_{folder}'
        if temp_folder.is_dir():
            shutil.rmtree(temp_folder)

    return


def warning_messange(task_name, run_folder, ymal_name):

    if 'arc' in task_name:
        with open(run_folder/f'{ymal_name}.yml', "r") as stream:
            yml_dict = yaml.safe_load(stream)
            for fits_frame in yml_dict['frames']:
                print('-- ', fits_frame, fits.getval(run_folder/f'data/{fits_frame}', 'VPH'), fits.getval(run_folder/f'data/{fits_frame}', 'SPECLAMP'))

    return


def indexing_the_frames(log, OB_task, VPH_task, type_task, obj_task=None, exclude_list=[]):

    if obj_task is None:
        idcs_frames = (log.OB == OB_task) & (log.VPH == VPH_task) & (log.type == type_task) & (~log.index.isin(exclude_list))
    else:
        idcs_frames = (log.OB == OB_task) & (log.VPH == VPH_task) & (log.object == obj_task) & (~log.index.isin(exclude_list))

    return idcs_frames


def store_task(task_ID_ref, VPH_ref, output_folder, DF_files, task_name, DF_tasks, idx_task, conf_task, idcs_files,
               arc_types):

    # Files list
    type_files = DF_files.loc[idcs_files].index.values

    # Remove non valid arc files
    if task_ID_ref == 'arc':
        task_files = []
        for file in type_files:
            lamp_type = fits.getval(f'{output_folder}/data/{file}', 'SPECLAMP')
            if lamp_type == arc_types:
                task_files.append(file)
    else:
        task_files = list(type_files)

    # Generate task yml
    yml_dict = {'id': task_ID_ref,
                'mode': megaradrp_modes[task_name],
                'instrument': 'MEGARA',
                'frames': task_files}
    yml_dict.update(conf_task)

    # Save yml to a text file
    dict_adress = f'{output_folder}/{task_ID_ref}.yml'
    with open(dict_adress, 'w') as f:
        yaml.dump(yml_dict, f, sort_keys=False)

    print(f'\n --- {task_ID_ref} ({np.sum(idcs_files)} files)--- ')
    for file in task_files:
        if task_ID_ref == 'arc':
            lamp_type = fits.getval(f'{output_folder}/data/{file}', 'SPECLAMP')
            message = f'{file} ({lamp_type})'
        else:
            message = f'{file}'
        print(message)

    # # Store DataFrame data
    # DF_tasks.loc[idx_task, 'task_id'] = task_ID_ref
    # DF_tasks.loc[idx_task, 'file_name'] = f'{task_ID_ref}.yml'
    # DF_tasks.loc[idx_task, 'VPH'] = VPH_ref

    return