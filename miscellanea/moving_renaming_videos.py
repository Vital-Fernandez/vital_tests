import os
import glob
import enzyme
from pathlib import Path
import shutil
import unidecode
import re
import string

videoFolder = Path('J:\Anime\Pocket Monsters') # TODO change to open format to avoid new dependency

forbidden_characters = "[(){}<>]',:?¿!¡"
file_names = os.listdir(videoFolder)
file_origNames = []
rootOutputFolder = Path('J:/Anime/myPokemon')

for video_file in file_names:
    with open(videoFolder / video_file, 'rb') as f:
        mkv = enzyme.MKV(f)
        video_data = mkv.info
        file_origNames.append(video_data.title)

move_check = True

for i, item in enumerate(file_names):

    # Declare original file
    err_name = file_names[i]
    inputFile = videoFolder / err_name
    file_name, extension = os.path.splitext(inputFile)

    # Metadata name
    orig_name = file_origNames[i]
    print(f'{i}: {inputFile}', os.path.isfile(inputFile))

    if orig_name is not None:

        fileReference = orig_name[8:orig_name.find('-')-1]
        seasonN = int(fileReference[1:fileReference.find('E')])
        episodeN = int(fileReference[fileReference.find('E')+1:])
        episodeTitle = orig_name[orig_name.find('-')+2:]
        seasonFolderName = r'Season {:02}'.format(seasonN)
        episodeName = r'S{:02}ES{:02} - {}{}'.format(seasonN, episodeN, episodeTitle, extension)

        #outputFolder = rootOutputFolder / seasonFolderName
        outputFolder = videoFolder / seasonFolderName
        outputFile = unidecode.unidecode(str(outputFolder / episodeName))
        outputFile = outputFile.translate(str.maketrans("", "", forbidden_characters))
        print('--- ', episodeName.translate(str.maketrans("", "", forbidden_characters)))

    #     # Create folder if not available
    #     Path(outputFolder).mkdir(parents=True, exist_ok=True)
    #
    #     # Move file if there is not conflict:
    #     if not os.path.isfile(outputFile):
    #         print(f'-- {outputFile}')
    #         if move_check:
    #             try:
    #                 # shutil.move(inputFile, outputFile)
    #                 os.rename(inputFile, outputFile)
    #                 print('-- Moved\n')
    #             except:
    #                 print('-- SHUTIL FAILURE')
    #     else:
    #         print(f'--  WARNING COULD NOT BE MOVED: {inputFile}\n')
    #
    # else:
    #     print(f'--  WARNING NO TITLE: {inputFile}\n')



