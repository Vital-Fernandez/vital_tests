import pathlib
import numpy as np
import pandas as pd
import src.specsiser as sr
import matplotlib.pyplot as plt
from astroquery.sdss import SDSS
import json
import wget
import time


def safe_dict(pet, filename):
    with open(filename, 'w') as f:
        f.write(json.dumps(pet))


def load_dict(filename):
    with open(filename) as f:
        pet = json.loads(f.read())
    return pet


def fetch_spec(out_file_address, plate, mjd, fiber):
    url = f'http://dr16.sdss.org/optical/spectrum/view/data/format=fits?plateid={plate}&mjd={mjd}&fiberid={fiber}&reduction2d=v5_7_0'
    wget.download(url, out_file_address)
    return


lineConversion = {'He2_4685A': 'He_II 4685',
                  'H1_4861A': 'H_beta',
                  'H1_6563A': 'H_alpha',
                  'O3_4363A': '[O_III] 4363',
                  'O3_4959A': '[O_III] 4959',
                  'O3_4959A': '[O_III] 5007',
                  'S3_6312A': '[S_III] 6312'}

conf_file_address = '../sampleHeII.ini'
obsData = sr.loadConfData(conf_file_address)

fits_folder = pathlib.Path(obsData['data_location']['fits_folder'])
data_folder = pathlib.Path(obsData['data_location']['treatment_folder'])
objList_file = obsData['data_location']['sample_file']
objDF = pd.read_csv(data_folder/objList_file, delim_whitespace=True, comment='#', names=['PLATEID', 'MDJ_1', 'FIBERID_1'])

normFlux = obsData['sample_data']['norm_flux']

# Loop through the sample list
err_dict = {}
n_objs = objDF.index.values.size

# Empty dataframe to store the data
# column_list = ['ID', 'plate', 'mjd', 'fiber', 'redshift', 'Hbeta', 'Hbeta_err', 'He2_4686A', 'He2_4686A_err']
# dataDF = pd.DataFrame(index=index, columns=)

# Prepare dataframe
objDF['index'] = objDF['PLATEID'].astype(str) + '-' + objDF['MDJ_1'].astype(str) + '-' + objDF['FIBERID_1'].astype(str)
objDF.set_index('index', inplace=True)

objDF['z'] = np.nan
for szLabel, sdssLabel in lineConversion.items():
    objDF[szLabel] = np.nan
    objDF[f'{szLabel}_err'] = np.nan
objDF['Comment'] = 'None'

for idx, obj in enumerate(objDF.index.values):

    # Object indexing
    plate, mjd, fiber = objDF.iloc[idx]['PLATEID':'FIBERID_1']
    objLabel = f'{plate}-{mjd}-{fiber}'
    spec_name = f'{objLabel}.fits'
    output_file = fits_folder / spec_name
    print(f'- Obj {idx}/{n_objs}')

    try:
        spec = SDSS.get_spectra(plate=plate, fiberID=fiber, mjd=mjd, timeout=60, cache=True, data_release=16)
        spec[0].writeto(output_file, overwrite=True)
    except:
        spec = None
        print(f'-- Astroquery failure: {objLabel}')

    if spec is None:
        try:
            output_file = fits_folder/spec_name
            fetch_spec(str(output_file), plate, mjd, fiber)
            time.sleep(3)
        except:
            spec = None
            objDF.loc[objLabel, 'Comment'] = 'No result query'
            err_dict[objLabel] = 'No result query'
            print(f'-- Manual download failure: {objLabel}')

        # Save the spectrum
    # if (spec is not None) and output_file.is_file():
    #     output_file = fits_folder/spec_name
    #     spec[0].writeto(output_file, overwrite=True)
    #
    #     # Get line data
    #     wave, data, hdrs = sr.import_fits_data(output_file, instrument='SDSS')
    #     dataProp = spec[0][2].data
    #     dataLines = spec[0][3].data
    #     objDF.loc[objLabel, 'z'] = dataProp['z'][0]
    #
    #     for szLabel, sdssLabel in lineConversion.items():
    #         idx_line = np.where(dataLines['LINENAME'] == sdssLabel)[0][0]
    #         objDF.loc[objLabel, szLabel] = dataLines['LINEAREA'][idx_line]
    #         objDF.loc[objLabel, f'{szLabel}_err'] = dataLines['LINEAREA_ERR'][idx_line]

    # Recover line data:
    # wave, data, hdrs = sr.import_fits_data(output_file, instrument='SDSS')
    # flux = data['flux'] * normFlux

# Print failing objects
print()
for key, value in err_dict.items():
    print(f'{key} \t {value}')

err_dict_file = data_folder/f'AVO_query_err_file.txt'
safe_dict(err_dict, err_dict_file)
myDict = load_dict(err_dict_file)

df_file = data_folder/f'AVO_dataframe.txt'
with open(df_file, 'wb') as output_file:
    string_DF = objDF.to_string()
    output_file.write(string_DF.encode('UTF-8'))

# with open(err_dict_file, 'w') as f:
#     # f.write(str(err_dict))
#     print(err_dict, file=f)

    # Plot the spectrum
    # wave, data, hdrs = sr.import_fits_data(output_file, instrument='SDSS')
    # flux = data['flux'] * normFlux
    # z_i = hdrs[1]["z"][0]
    # lm = sr.LineMesurer(wave, flux, redshift=z_i, normFlux=normFlux)
    # lm.plot_spectrum_components()