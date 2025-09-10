import os
import numpy as np
from pathlib import Path

file_list = Path(r'D:\Cofre\Photos\files.txt')
root_address = Path(r'D:\Cofre\Photos')

with open(file_list) as f:
    photo_list = f.readlines()

photo_list = np.array(photo_list)

for photo_address in photo_list:
    photo_path = root_address/photo_address[:-1]
    if ('Orb Library' in photo_address) and (photo_path.is_file()):
        print(photo_path, photo_path.is_file())
        # photo_path.unlink()
