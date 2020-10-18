import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path

objList = ['gp030321', 'gp101157', 'gp121903']
conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList=objList, group_variables=False)
starlightFolder = Path(obsData['SSP_synthesis']['starlight_folder'])

fileList = obsData['file_information']['files_list']
dataFolder = Path(obsData['file_information']['data_folder'])
outputFolder = dataFolder/'flux_analysis'

objList_B = obsData['file_information']['objectB_list']
fileList_B = obsData['file_information']['filesB_list']
objList_R = obsData['file_information']['objectR_list']
fileList_R = obsData['file_information']['filesR_list']

z_objs = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
flux_norm = obsData['sample_data']['norm_flux']
noise_region = obsData['sample_data']['noiseRegion_array']
idx_band = int(obsData['file_information']['band_flux'])

counter = 0
for i, obj in enumerate(objList):

    z = z_objs[i]
    cHbeta = obsData[obj]['cHbeta']

    for ext in ('_BR', '_B'):

        # # Establish files location
        # objName = obsData['file_information']['object_list'][i]
        # fitsFolder, fitsFile = file_address.parent, file_address.name
        # masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')
        # lineLogFolder, lineLogFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_linesLog_Emission.txt')
        # emissionSpectrumAddress = fitsFolder/'flux_analysis'/f'{objName}_gasSpectrum.txt'
        # plotFolder, plotFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_linesGrid_Emission')

        # Declare files location
        fits_file = dataFolder/f'{obj}{ext}.fits'
        lineLog_file = outputFolder/f'{obj}{ext}_linesLog_Emission.txt'
        objMask = dataFolder/'flux_analysis'/f'{obj}{ext}_mask.txt'
        objGasSpectrumFile = outputFolder/f'{obj}{ext}_gasSpectrum.txt'
        fit_conf = obsData[f'{obj}_line_fitting']
        lineGrid_file = dataFolder/'flux_analysis'/f'{obj}{ext}_linesGrid_Emission.png'

        # Get fits data
        wave, flux = np.loadtxt(objGasSpectrumFile, unpack=True)

        # Load line measurer object
        maskDF = pd.read_csv(objMask, delim_whitespace=True, header=0, index_col=0)
        lm = sr.LineMesurer(wave, flux, objMask, normFlux=flux_norm)

        # Loop through the lines
        print(f'\n-- Treating {counter} :{obj}{ext}.fits')

        # Fit and check the regions
        obsLines = maskDF.index.values
        for j, lineLabel in enumerate(obsLines):

            print(f'-- {lineLabel}:')
            wave_regions = maskDF.loc[lineLabel, 'w1':'w6'].values
            lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf=fit_conf)

            if lm.blended_check:
                plotFile = f'{obj}{ext}_{lineLabel}_emission.png'
                lm.plot_fit_components(lm.fit_output, output_address=outputFolder/plotFile)
        lm.save_lineslog(lm.linesDF, lineLog_file)

        # Plot the single lines:
        idcs_unblended = ~lm.linesDF.index.str.contains('_b')
        lm.plot_line_grid(lm.linesDF.loc[idcs_unblended], ncols=8, output_address=lineGrid_file)

        # Increase counter for obj number
        counter += 1
