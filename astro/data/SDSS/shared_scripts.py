import os


def list_objName(directory, extension):
    output_list = []
    for file in os.listdir(directory):
        if file.endswith(extension):
            output_list.append(os.path.splitext(file)[0])
    return output_list


obsFolder = 'D:/Google drive/Astrophysics/Datos/SDSS-Ricardo/green_peas/'
obsConfaddress = 'D:/Pycharm Projects/vital_tests/astro/data/SDSS/flux_comparison.ini'
