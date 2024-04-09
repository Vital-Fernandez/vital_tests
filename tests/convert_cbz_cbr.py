import os
import shutil
import unrar
from unrar import rarfile
import zipfile


def convert_rar_to_zip(folder_path):
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)
            if file.endswith(".cbr"):
                try:
                    with rarfile.RarFile(file_path) as rf:
                        rf.extractall(root)
                        for extracted_file in rf.namelist():
                            extracted_file_path = os.path.join(root, extracted_file)
                            if os.path.isfile(extracted_file_path):
                                zip_file_path = os.path.splitext(extracted_file_path)[0] + '.cbz'
                                with zipfile.ZipFile(zip_file_path, 'w', zipfile.ZIP_DEFLATED) as zf:
                                    zf.write(extracted_file_path, os.path.basename(extracted_file_path))
                                os.remove(extracted_file_path)
                except Exception as e:
                    print(f"Error processing {file_path}: {str(e)}")


if __name__ == "__main__":
    folder_path = "/home/vital/Desktop"  # Replace with the path to your folder
    convert_rar_to_zip(folder_path)