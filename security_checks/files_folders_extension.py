import pathlib
str_fileAddres = "/Users/pankaj/abc.txt"
fileAddres = pathlib.Path(str_fileAddres)
fileAddres2 = pathlib.Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
print(fileAddres.suffix)
print(fileAddres.name)
print(fileAddres.parent)
print(fileAddres.parts)
print(pathlib.Path(str_fileAddres).with_suffix('.db'))
print(fileAddres2.is_file())
print(fileAddres2.parent.exists())