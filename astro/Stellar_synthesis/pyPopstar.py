import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import zipfile
from pathlib import Path

SSP_zip = Path(f'S:\Astro_data\Observations\HR-pypostar\CHA_Z004.zip')

# root_zip = zipfile.Path(models_path)
#
# print(list(root_zip.iterdir()))
#
# directory = pathlib.Path("source_dir/")
#
# with zipfile.ZipFile(SSP_zip, 'r') as zip_ref:
#     for file in zip_ref.namelist():
#         print(file)
print(SSP_zip.is_file())
with zipfile.ZipFile(SSP_zip.as_posix(), 'r') as zip_ref:
    directories = set()
    for zip_info in zip_ref.infolist():
        if zip_info.filename.endswith('/'):  # check if it is a directory
            directories.add(zip_info.filename)

   # for file_path in directory.iterdir():
   #     archive.write(file_path, arcname=file_path.name)