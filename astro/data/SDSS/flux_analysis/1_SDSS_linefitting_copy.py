import os
import numpy as np
import pandas as pd
import src.specsiser as sr
from matplotlib import pyplot as plt, rcParams
from pathlib import Path


def list_files(directory, extension):
    output_list = []
    for file in os.listdir(directory):
        if file.endswith(extension):
            output_list.append(os.path.join(directory, file))
    return output_list


# Import the observation data
obsData = sr.loadConfData('D:/Pycharm Projects/vital_tests/astro/data/SDSS/flux_comparison.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
linesDb = pd.read_excel(linesFile, sheet_name=0, header=0, index_col=0)
data_folder = Path(obsData['file_information']['data_folder'])
fileList = list_files(data_folder, '.fits')
addressList = list(data_folder / file for file in fileList)
flux_norm = obsData['sample_data']['norm_flux']

# # Analyse the spectrum
for i, file_address in enumerate(addressList):

    if i==1:

        # Open lineslog
        fitsFolder, fitsFile = file_address.parent, file_address.name
        masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')
        lineLogFolder, lineLogFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_linesLog.txt')
        plotFolder, plotFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_singleLines')
        objName = fitsFile.replace('.fits', '')

        # Set and crop the wavelength
        wave_rest, flux, header = sr.import_fits_data(fitsFolder/fitsFile, instrument='SDSS')
        idx_wave = (wave_rest >= obsData['sample_data']['wmin_array']) & (wave_rest <= obsData['sample_data']['wmax_array'])

        # Load line measurer object
        lm = sr.LineMesurerGUI(wave_rest[idx_wave], flux[idx_wave], masksFolder/masksFile, normFlux=flux_norm)

        # Loop through the lines
        print(f'\n- {i}: {objName}')
        obsLines = lm.linesDF.index.values
        for j, lineLabel in enumerate(obsLines):

            if lineLabel == 'H1_6563A_b':

                # Declare regions data
                print(f'-- {lineLabel}:')
                wave_regions = lm.linesDF.loc[lineLabel, 'w1':'w6'].values
                idcsLinePeak, idcsContinua = lm.define_masks(wave_regions)

                # Measure line data
                lm.line_properties(idcsLinePeak, idcsContinua, bootstrap_size=1000)

                # Establish line and object fit configuration # TODO put all this tricks in an external method
                fit_conf = {}
                if lineLabel in obsData['blended_groups']:
                    fit_conf[lineLabel] = obsData['blended_groups'][lineLabel]
                if 'default_blended_lines' in obsData:
                    fit_conf.update(obsData['default_blended_lines'])
                if f'{objName}_blended_lines' in obsData:
                    fit_conf.update(obsData[f'{objName}_blended_lines'])

                for key, value in fit_conf.items():
                    print(key, value)

                fig, ax = plt.subplots()
                ax.step(lm.wave[idcsLinePeak], lm.flux[idcsLinePeak], label='Line spectrum')
                ax.step(lm.wave[idcsContinua], lm.flux[idcsContinua], label='Line continuum')
                ax.legend()
                plt.show()

                # Fit the emission lines
                lm.line_fit('lmfit', lineLabel, idcsLinePeak, idcsContinua, continuum_check=True, user_conf=fit_conf)

                # Save deblending fitting
                if lineLabel in obsData['blended_groups']:
                    plotBlendedFile = fitsFile.replace('.fits', f'_{lineLabel}_deblending')
                    lm.plot_fit_components(lm.fit_output)



    # # Save dataframe to text file
    # lm.linesDF.sort_values('mu', inplace=True)
    # lm.save_lineslog(lm.linesDF, lineLogFolder/lineLogFile)
    #
    # # Plot the single lines:
    # lm.plot_detected_lines(lm.linesDF, ncols=5, output_address=plotFolder/plotFile)