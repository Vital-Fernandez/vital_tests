import os
import glob

downloads_folder = '/home/vital/Downloads/Astrodata/20110925/'
output_folder = '/home/vital/Astrodata/WHT_2011_09/Night2/raw_fits/'

print(os.listdir(downloads_folder))
print(glob.glob(downloads_folder))
print(glob.glob('/home/vital/Downloads/Astrodata/'))
print(glob.glob(downloads_folder+ '*.fits.fz'))

list_files = map(( lambda x: downloads_folder + x), os.listdir(downloads_folder))

for fz_file in os.listdir(downloads_folder):
    input_file      = downloads_folder + fz_file
    extract_command = 'funpack {}'.format(input_file)
    os.system(extract_command)

for fz_file in os.listdir(downloads_folder):
    uncompress_file = fz_file.replace('.fz', '')
    output_name     = 'r' + uncompress_file[uncompress_file.find('_')+1:len(uncompress_file)][1:]
    if os.path.isfile(downloads_folder + uncompress_file):
        os.rename(downloads_folder + uncompress_file, output_folder + output_name)
    else:
        print 'This file fails', uncompress_file