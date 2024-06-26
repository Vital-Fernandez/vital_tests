import numpy as np
import lime as lm
from astropy.io import fits
from pathlib import Path
from astro.papers.SHOC579_project.MEGARA_reduction.INAOE_tools.convert import convert

obs_conf = lm.load_cfg(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\obsConf.ini')
z_obj = obs_conf['sample_data']['z_array'][0]
norm_flux = obs_conf['sample_data']['norm_flux']

reduced_folder = Path('S:/Astro_data/Observations/SHOC579/MEGARA/obsidOB0001_LR-B_5_SHOC579_lcb_image_6_result/')

final_spec_address = reduced_folder/'final_rss.fits'
reduced_spec_address = reduced_folder/'reduced_image.fits'
convert_spect_adddress = reduced_folder/'convert_rss.fits'

# fits.info(reduced_spec_address)
# fits.info(final_spec_address)

# hdu_list = fits.open(reduced_spec_address)
# data0, header0 = hdu_list[0].data, hdu_list[0].header
# data1, header1 = hdu_list[1].data, hdu_list[1].header
# data2, header2 = hdu_list[2].data, hdu_list[2].header
# hdu_list.close()

# Tool from Javier
cube_hdu = convert(final_spec_address, writeout=convert_spect_adddress, overwrite=True)

# plot_spectrum()

cube_hdu = fits.open(convert_spect_adddress)
data = cube_hdu[0].data
header = cube_hdu[0].header

dw = header['CD3_3']
w_min = header['CRVAL3']
nPixels = header['NAXIS3']
w_max = w_min + dw * nPixels
wave = np.linspace(w_min, w_max, nPixels, endpoint=False)

# Hbeta map
line_region = np.array([5080, 5100])
idcsLine = np.searchsorted(wave, line_region)
line_map = data[idcsLine[0]:idcsLine[1], :, :].sum(axis=0)


lm.CubeFitsInspector(wave, data, image_bg=line_map, norm_flux=norm_flux, redshift=z_obj)

# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot()
# im = ax.imshow(line_map, interpolation='none')
# ax.update({'title': r'SHOC579 galaxy', 'xlabel': r'X', 'ylabel': r'Y'})
# plt.show()
#
#
# idx_j, idx_i = 20, 22
#
# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot()
# im = ax.step(wave, data[:, idx_j, idx_i])
# ax.update({'title': r'SHOC579 galaxy', 'xlabel': r'Wave', 'ylabel': r'flux'})
# plt.show()


# visualization.main([str(final_spec_address)])

# fig = plt.figure(figsize=(10, 10))
# ax = fig.add_subplot()
# im = ax.imshow(np.log10(data0))
# ax.update({'title': r'J0926 galaxy', 'xlabel': r'X', 'ylabel': r'Y'})
# plt.show()




# # Loading project configuration file
# obs_conf = lm.load_cfg(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\obsConf.ini')
# reduction_cfg = obs_conf['Megara_reduction']
# obj_list = reduction_cfg['obj_list']
# std_list = reduction_cfg['std_star_list']
#
# # Dataframe with files list
# rd_df_address = Path(reduction_cfg['rd_df_address'])
# sample_DF = lm.load_lines_log(f'{rd_df_address}.txt')
#
# # Stating folder structure
# instructions_folder = Path(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\MEGARA_reduction')
# reduction_folder = Path(reduction_cfg['root_folder'])
# data_folder = reduction_folder/'data'
#
# # Run the pipeline one OB at a time:
# OB_list = ['OB0001']
# complete_task_list = ['bias', 'trace_map', 'arc', 'fiber_flat', 'lcb_std', 'lcb_image']
#
# for OB in OB_list:
#
#     # Create clean requirements file at each run
#     original_yml = instructions_folder/'SHOC579_req.yml'
#     req_yml = reduction_folder/f'control_{OB}.yaml'
#     shutil.copyfile(original_yml, req_yml)
#
#     # Define data manager
#     dm = create_datamanager(req_yml, reduction_folder, data_folder)
#
#     # Generate complete list of tasks
#     tasksNames_dict = {}
#     for idx_task, task in enumerate(['bias', 'trace_map', 'arc', 'fiber_flat']):
#         tasksNames_dict[task] = f'{OB}_{task}'
#     for std in std_list:
#         for task in ['lcb_std', 'lcb_image']:
#             if Path(reduction_folder/f'{OB}_{std}_{task}.yml').is_file():
#                 tasksNames_dict[f'{std}_{task}'] = f'{OB}_{std}_{task}'
#     for obj in obj_list:
#         task = 'lcb_image'
#         if Path(reduction_folder / f'{OB}_{obj}_{task}.yml').is_file():
#             tasksNames_dict[f'{obj}_{task}'] = f'{OB}_{obj}_{task}'
#
#     # delete_task_temp_folder(OB, task, idx_task, reduction_folder)
#
#     # Load the observation files
#     OB_task_list = list(tasksNames_dict.values())
#     task_file_list = [task + '.yml' for task in OB_task_list]
#     with ctx.working_directory(reduction_folder):
#         sessions, loaded_obs = load_observations(task_file_list, is_session=False)
#         dm.backend.add_obs(loaded_obs)
#
#     # Run the tasks
#     task_list = list(tasksNames_dict.values())
#     for run_id in task_list:
#         print(f'Running: {run_id}')
#         output_run = run_reduce(dm, run_id)
