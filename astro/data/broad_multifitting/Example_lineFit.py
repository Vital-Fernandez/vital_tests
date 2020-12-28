from pathlib import Path
import src.specsiser as sr
from fitelp.read_spectra import read_spectra

data_folder = Path('D:/Google drive/Astrophysics/Datos/broad_multiComponent')
obsData = sr.loadConfData('broad_conf.ini', group_variables=False)
z_mean, z_err = sr.redshift_calculation(obsData['sample_data']['obs_waves'], obsData['sample_data']['emis_waves'])
norm_flux = obsData['sample_data']['norm_flux']
obj_list = ['B6479s', 'R8731s']

for idx_obj, obj in enumerate(['R8731s']):

    fits_address = data_folder/f'{obj}.fits'
    wave_data, flux_data = read_spectra(fits_address, scaleFlux=1)

    mask_address = data_folder/f'{obj}_mask.txt'
    mask_DF = sr.lineslogFile_to_DF(mask_address)

    # Individual line measurement
    lm = sr.LineMesurer(wave_data[0],  flux_data[0], normFlux=norm_flux, redshift=z_mean)
    fitConf = obsData[f'default_line_fitting']

    # Loop through the lines
    for lineLabel in ['H1_6563A_b']:

        # Get line ranges
        lineWaves = mask_DF.loc[lineLabel, 'w1':'w6'].values

        # Perform fit
        lm.fit_from_wavelengths(lineLabel, lineWaves, fit_conf=fitConf)

        # Display results
        lm.print_results(show_fit_report=True, show_plot=True)


