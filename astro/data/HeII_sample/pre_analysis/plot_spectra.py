import pathlib
import numpy as np
import pandas as pd
import src.specsiser as sr

conf_file_address = '../sampleHeII.ini'
obsData = sr.loadConfData(conf_file_address)

fits_folder = pathlib.Path(obsData['data_location']['fits_folder'])
data_folder = pathlib.Path(obsData['data_location']['treatment_folder'])
objList_file = data_folder/f'AVO_dataframe.txt'

pertil_array = obsData['sample_data']['percentil_array']

# Load the list of objects
logDF = sr.lineslogFile_to_DF(objList_file)

arrayFluxLevel = np.percentile(logDF['He2_4685A'].values, pertil_array)

# Recover the group of objects
for i, flux_level in enumerate(arrayFluxLevel):

    if i == 0:
        idcs_spec = logDF['He2_4685A'] >= flux_level
    else:
        idcs_spec = (logDF['He2_4685A'] < arrayFluxLevel[i-1]) & (logDF['He2_4685A'] >= flux_level)

    logDF.loc[idcs_spec, 'intensity_Group'] = i

df_file = data_folder/f'AVO_catalogue_dataframe.txt'
with open(df_file, 'wb') as output_file:
    string_DF = logDF.to_string()
    output_file.write(string_DF.encode('UTF-8'))

print(logDF)

    # # Recover the spectrum
    # wave, data, hdrs = sr.import_fits_data(sampleFiles[i], instrument='SDSS')
    #
    # print(i, flux_level, np.sum(idcs_spec))
