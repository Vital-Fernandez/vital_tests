import src.specsiser as sr
from pathlib import Path
from astro.data.SDSS.shared_scripts import list_objName, obsConfaddress, obsFolder


# Declare data and files location
objList = list_objName(obsFolder, '.fits')
obsData = sr.loadConfData(obsConfaddress, objList, group_variables=False)
data_folder = Path(obsData['file_information']['data_folder'])
addressList = list(Path(f'{data_folder/objName}.fits') for objName in objList)

# Sample properties
norm_Flux = obsData['sample_data']['norm_flux']
wmin, wmax = obsData['sample_data']['wmin_array'], obsData['sample_data']['wmax_array']

# Analyse the spectrum
for i, file_address in enumerate(addressList):

    # if i == 0:

        # Open lineslog
        objName = objList[i]
        fitsFolder, fitsFile = file_address.parent, file_address.name
        masksFolder, masksFile = fitsFolder, fitsFile.replace('.fits', '_masks.txt')
        lineLogFolder, lineLogFile = fitsFolder/'flux_analysis', fitsFile.replace('.fits', '_linesLog.txt')
        plotFolder, plotGrid = fitsFolder / 'flux_analysis', fitsFile.replace('.fits', '_lineGrid.png')
        print(f'\n- {i}: {objName}')

        # Set and crop the wavelength
        wave_rest, flux, header = sr.import_fits_data(fitsFolder/fitsFile, instrument='SDSS')
        idx_wave = (wave_rest >= wmin) & (wave_rest <= wmax)

        # Load line measurer object
        lm = sr.LineMesurer(wave_rest[idx_wave], flux[idx_wave], masksFolder/masksFile, normFlux=norm_Flux)

        # Get lines and their fitting configuration
        obsLines = lm.linesDF.index.values
        fit_conf = obsData[f'{objName}_line_fitting']

        # Loop through the lines
        for j, lineLabel in enumerate(obsLines):

            # if 'H1_6563' in lineLabel:

            # Fit each line regions data
            print(f'-- {lineLabel}:')
            wave_regions = lm.linesDF.loc[lineLabel, 'w1':'w6'].values
            lm.fit_from_wavelengths(lineLabel, wave_regions, fit_conf)

            if lm.blended_check:
                plotBlend = fitsFile.replace('.fits', f'{lineLabel}_deblending.png')
                lm.plot_fit_components(lm.fit_output, output_address=plotFolder/plotBlend)

        lm.save_lineslog(lm.linesDF, lineLogFolder/lineLogFile)

        # Plot the single lines:
        idcs_unblended = ~lm.linesDF.index.str.contains('_b')
        lm.plot_line_grid(lm.linesDF.loc[idcs_unblended], ncols=8, output_address=plotFolder/plotGrid)

