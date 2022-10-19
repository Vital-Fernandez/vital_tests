import os
import pathlib
import logging
import tarfile

import numina.util.context as ctx
from numina.user.helpers import create_datamanager, load_observations
from numina.user.baserun import run_reduce
from numina.tests.testcache import download_cache


import matplotlib.pyplot as plt
logging.basicConfig(level=logging.DEBUG)
basedir = pathlib.Path().resolve()
tarball = 'MEGARA-cookbook-M15_LCB_HR-R-v1.tar.gz'
url = 'http://guaix.fis.ucm.es/~spr/megara_test/{}'.format(tarball)

downloaded = download_cache(url)

# Uncompress
with tarfile.open(downloaded.name, mode="r:gz") as tar:
    def is_within_directory(directory, target):
        
        abs_directory = os.path.abspath(directory)
        abs_target = os.path.abspath(target)
    
        prefix = os.path.commonprefix([abs_directory, abs_target])
        
        return prefix == abs_directory
    
    def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
    
        for member in tar.getmembers():
            member_path = os.path.join(path, member.name)
            if not is_within_directory(path, member_path):
                raise Exception("Attempted Path Traversal in Tar File")
    
        tar.extractall(path, members, numeric_owner=numeric_owner) 
        
    
    safe_extract(tar)


# import pathlib
# import yaml
# import shutil
#
# from numina.util import context as ctx
# from numina.user.helpers import create_datamanager, load_observations
# from numina.user.baserun import run_reduce
#
# # Stating folder structure
# instructions_folder = pathlib.Path(r'D:\Pycharm Projects\vital_tests\astro\papers\SHOC579_project\MEGARA_reduction')
# reduction_folder = pathlib.Path(r'S:\Astro_data\Observations\MEGARA_test')
# data_folder = reduction_folder/'data'
#
# # Paste clean requirements file at each run
# original_yml = instructions_folder/'control_orig.yaml'
# req_yml = reduction_folder/'control_v2.yaml'
# shutil.copyfile(original_yml, req_yml)
#
# # Define data manager
# dm = create_datamanager(req_yml, reduction_folder, data_folder)
#
# # Load observational files at reduction folder
# with ctx.working_directory(reduction_folder):
#     obsresults = ["0_bias.yaml",
#                   "1_M15_tracemap.yaml",
#                   "2_M15_modelmap.yaml",
#                   "3_M15_wavecalib.yaml",
#                   "4_M15_fiberflat.yaml",
#                   "5_M15_twilight.yaml",
#                   "6_M15_Lcbadquisition.yaml",
#                   "8_M15_reduce_LCB.yaml",
#                   "7_M15_Standardstar.yaml"]
#
#     sessions, loaded_obs = load_observations(obsresults, is_session=False)
#     dm.backend.add_obs(loaded_obs)
#
# # Run the tasks
# obsid = "0_bias"
# print(f'- Running {obsid}')
# task0 = run_reduce(dm, obsid)
#
# obsid = "1_HR-R"
# print(f'- Running {obsid}')
# task1 = run_reduce(dm, obsid)
#
# obsid = "2_HR-R"
# print(f'- Running {obsid}')
# task1 = run_reduce(dm, obsid)

# obsid = "3_HR-R"
# print(f'- Running {obsid}')
# task3 = run_reduce(dm, obsid)

# obsid = "4_HR-R"
# print(f'- Running {obsid}')
# task4 = run_reduce(dm, obsid)
#
# obsid = "5_HR-R"
# print(f'- Running {obsid}')
# task5 = run_reduce(dm, obsid)
#
# obsid = "6_HR-R"
# print(f'- Running {obsid}')
# task6 = run_reduce(dm, obsid)
#
# obsid = "7_HR-R"
# print(f'- Running {obsid}')
# task7 = run_reduce(dm, obsid)
#
# obsid = "8_HR-R"
# print(f'- Running {obsid}')
# task8 = run_reduce(dm, obsid)

# Back up of task-updated requirement file
with ctx.working_directory(reduction_folder):
    with open('control_dump.yaml', 'w') as fd:
        datam = dm.backend.dump_data()
        yaml.dump(datam, fd)

# id: 4_HR-R
# mode: MegaraFiberFlatImage
#
# id: 5_HR-R
# mode: MegaraTwilightFlatImage
#
# id: 6_HR-R
# mode: MegaraLcbAcquisition
#
# id: 7_HR-R
# mode: MegaraLcbStdStar
#
# id: 8_HR-R
# mode: MegaraLcbImage
#
