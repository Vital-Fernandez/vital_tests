from pathlib import Path
import numpy as np
import src.specsiser as sr

# Import the observation data
obsData = sr.loadConfData('../gtc_greenpeas_data.ini', group_variables=False)
linesFile = Path('D:/Pycharm Projects/spectra-synthesizer/src/specsiser/literature_data/lines_data.xlsx')
data_folder = Path(obsData['file_information']['data_folder'])
file_list = obsData['file_information']['files_list']
addressList = list(data_folder/file for file in file_list)
flux_norm = obsData['sample_data']['norm_flux']

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # Establish files location
    objName = obsData['file_information']['object_list'][i]
    fitsFolder, fitsFile = file_address.parent, file_address.name
    masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')
    lineLogFolder, lineLogFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_linesLog_Emission.txt')
    emissionSpectrumAddress = fitsFolder/'flux_analysis'/f'{objName}_gasSpectrum.txt'
    plotFolder, plotFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_linesGrid_Emission')

    # Get fits data
    wave, flux = np.loadtxt(emissionSpectrumAddress, unpack=True)

    # Load line measurer object
    lm = sr.LineMesurerGUI(wave, flux, masksFolder/masksFile, normFlux=flux_norm)

    # Loop through the lines
    print(f'\n- {objName}:')
    obsLines = lm.linesDF.index.values
    for j, lineLabel in enumerate(obsLines):

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
        if f'{objName}_blended_lines' in obsData:
            fit_conf.update(obsData[f'{objName}_blended_lines'])

        # Fit the emission lines
        lm.line_fit('lmfit', lineLabel, idcsLinePeak, idcsContinua, continuum_check=True, user_conf=fit_conf)

        # Save deblending fitting
        if lineLabel in obsData['blended_groups']:
            plotBlendedFile = fitsFile.replace('.fits', f'_{lineLabel}_deblending')
            # lm.plot_fit_components(lm.fit_output, output_address=plotFolder/plotBlendedFile)

    # Save dataframe to text file
    lm.linesDF.sort_values('wavelength', inplace=True)
    lm.save_lineslog(lm.linesDF, lineLogFolder/lineLogFile)

    # Plot the single lines:
    idcs_unblended = ~lm.linesDF.index.str.contains('_b')
    lm.plot_detected_lines(lm.linesDF.loc[idcs_unblended], ncols=8, output_address=plotFolder/plotFile)
