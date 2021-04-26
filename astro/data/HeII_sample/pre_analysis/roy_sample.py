import pathlib
import src.specsiser as sr
import numpy as np
import lineid_plot
import matplotlib.pyplot as plt

conf_file_address = '../sampleHeII.ini'
obsData = sr.loadConfData(conf_file_address)

dataFolder = pathlib.Path(obsData['data_location']['data_folder'])
treatmentFolder = pathlib.Path(obsData['data_location']['treatment_folder'])
resultsFolder = pathlib.Path(obsData['data_location']['results_folder'])

sampleFiles = tuple(dataFolder.iterdir())
sampleObj = tuple([x.name[x.name.find('-')+1:x.name.find('-')+8] for x in sampleFiles])
normFlux = obsData['sample_data']['norm_flux']

for i, obj in enumerate(sampleObj):

    # Load the data
    wave, data, hdrs = sr.import_fits_data(sampleFiles[i], instrument='SDSS')
    flux = data['flux'] * normFlux
    z_i = hdrs[1]["z"][0]

    #Output files
    specPlotAddress = f'{obj}_SDSS_spectrum.png'

    print(f'\nGalaxy {obj}')
    lineLabels = ['He_II 4685', 'H_beta', 'H_alpha', '[O_III] 4959', '[O_III] 5007']
    lineWaves = [4685.0, 4861.0, 5007]
    for lineRef in lineLabels:
        idx_line = np.where(hdrs[2]['LINENAME'] == lineRef)[0][0]
        lineArea, lineAreaErr = hdrs[2]['LINEAREA'][idx_line], hdrs[2]['LINEAREA_ERR'][idx_line]
        print(f'{lineRef} linea area : {lineArea:.2f}+/-{lineAreaErr:.2f}')

    lm = sr.LineMesurer(wave, flux, redshift=z_i, normFlux=normFlux, crop_waves=(1+z_i) * np.array([4685-100, 5100]))
    # lm.plot_spectrum_components(specLabel=f'Galaxy {obj}',
    #                             axConf={'ylabel': r'Flux $(10^{17}\,erg\,cm^{-2} s^{-1} \AA^{-1})$'},
    #                             output_address=resultsFolder/specPlotAddress)
    # lm.plot_spectrum_components(specLabel=f'Galaxy {obj}',
    #                             axConf={'ylabel': r'Flux $(10^{17}\,erg\,cm^{-2} s^{-1} \AA^{-1})$'})


    lineWaves = [4685.0, 4861.0, 5007]
    lineLabels = [r'He_II4685\AA', r'$H\beta$', r'[O_III]5007\AA']

    lineid_plot.plot_line_ids(lm.wave, lm.flux, lineWaves, lineLabels)

    plt.show()