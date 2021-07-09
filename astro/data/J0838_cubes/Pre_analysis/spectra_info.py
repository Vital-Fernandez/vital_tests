from pathlib import Path
import src.specsiser as sr
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt, rcParams, gridspec
from src.specsiser.tools.line_fitting import EmissionFitting
from astro.data.muse.common_methods import compute_line_flux_image, lineAreas, IFU_Cube_Plotter
from astro.data.J0838_cubes.common_methods import import_fits_data

# Declare data and files location
obsConf = sr.loadConfData('../J0838_cubes.ini')
fitsFolder = Path(obsConf['data_location']['fits_folder'])
dataFolder = Path(obsConf['data_location']['data_folder'])
resultsFolder = Path(obsConf['data_location']['results_folder'])

fileList = obsConf['data_location']['file_list']
objList = obsConf['data_location']['object_list']
z_list = obsConf['sample_data']['z_array']
norm_flux = obsConf['sample_data']['norm_flux']
ef = EmissionFitting()
percentil_array = obsConf['sample_data']['percentil_array']

for i, obj in enumerate(objList):

    if i == 0:

        # Input data
        objFolder = resultsFolder
        cube_address = fitsFolder / fileList[i]
        mask_global_DF = sr.lineslogFile_to_DF(dataFolder/obsConf['data_location']['mask_global'])

        # Output data
        db_addresss = objFolder / f'{obj}_database.fits'
        hdul_lineslog = fits.HDUList()

        # Load the data
        print('\n', db_addresss)
        wave, data, header = import_fits_data(cube_address, crval3=obsConf[obj]['CRVAL3'], frame_idx=0)

        for lineLabel in ['O3_5007A_b']:
            lineWaves = mask_global_DF.loc[lineLabel, 'w1':'w6']
            idcsLineRegion, idcsContRegion = ef.define_masks(wave/(1+z_list[i]), data, lineWaves)
            line_slice = data[idcsLineRegion, :, :]
            line_cont = data[idcsContRegion, :, :]
            avg = np.sum(line_slice, axis=0)

            rms = np.sqrt(np.mean(np.power(line_slice - avg, 2), axis=0))
            SNR_image = avg / rms

            flux_levels = np.percentile(SNR_image, percentil_array)

            IFU_Cube_Plotter(wave, data, line_slice.sum(axis=0),  header=header)


            # fig = plt.figure(figsize=(18, 5))
            # sky_wcs = WCS(header)
            # # ax = fig.add_subplot(projection=sky_wcs, slices=('x', 'y', 1))
            # ax = fig.add_subplot()
            # im = ax.imshow(line_cont.sum(axis=0))
            # plt.show()


        # -------------------------- Line mask generation --------------------------
        # # Voxel treatment
        # print(z_list[i])
        # flux_voxel = data[:, 8, 10]
        # lm = sr.LineMesurer(wave, flux_voxel, redshift=z_list[i], normFlux=norm_flux)
        # lm.plot_spectrum()
        #
        # # Find lines
        # norm_spec = lm.continuum_remover(noiseRegionLims=obsConf[obj]['noiseRegion_array'])
        # obsLinesTable = lm.line_finder(norm_spec, noiseWaveLim=obsConf[obj]['noiseRegion_array'], intLineThreshold=2.5)
        # matchedDF = lm.match_lines(obsLinesTable, mask_global_DF, tol=10, find_line_borders=False)
        # lm.plot_spectrum(obsLinesTable=obsLinesTable, matchedLinesDF=matchedDF, specLabel=f'Emission line detection')
        #
        # # Correct line region
        # # local_mask_DF = sr.lineslogFile_to_DF(dataFolder/'J0838_global_mask_red.txt')
        #
        # corrected_mask_file = Path(dataFolder/'B_mask_corrected.txt')
        # # lm.plot_line_mask_selection(matchedDF, corrected_mask_file, logscale=True)
        # lm.plot_line_grid(matchedDF, frame='obs', log_scale=False)
