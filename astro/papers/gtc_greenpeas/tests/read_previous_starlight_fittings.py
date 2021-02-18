import numpy as np
import os
from pathlib import Path
import src.specsiser as sr
from src.specsiser.physical_model.starContinuum_functions import SSPsynthesizer, computeSSP_galaxy_mass
from scipy.interpolate import interp1d
from astro.papers.gtc_greenpeas.common_methods import double_arm_redCorr
import pyneb as pn

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

objList = obsData['file_information']['object_list']
dataFolder = Path(obsData['file_information']['data_folder'])
resultsFolder = Path(obsData['file_information']['results_folder'])
fileList = obsData['file_information']['files_list']

fileList = obsData['file_information']['files_list']
idx_band = int(obsData['file_information']['band_flux'])

z_array = obsData['sample_data']['z_array']
wmin_array = obsData['sample_data']['wmin_array']
wmax_array = obsData['sample_data']['wmax_array']
arm_wave_boundary = obsData['sample_data']['w_div']

w_div_array = obsData['sample_data']['w_div']
red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

#This is the one
ext = 'BR'
cycle = 'Amorin2012'

for i, obj in enumerate(objList):

    if i > 2:

        print(f'\n-- Treating: {obj}{ext}.fits')

        # Declare input files
        objFolder = resultsFolder / f'{obj}'
        fits_file = dataFolder / f'{obj}_{ext}.fits'
        results_file = objFolder / f'{obj}_{ext}_measurements.txt'
        starlight2012_folder = dataFolder/'Starlight2012'
        outputFile = f'{obj}.BN'

        # Declare output files
        massFracPlotFile = objFolder / f'{obj}_{ext}_SSP_MasFrac_Papaderos_{cycle}.png'
        LightFracPlotFile = objFolder / f'{obj}_{ext}_SSP_LightFrac_Papaderos_{cycle}.png'
        stellarPlotFile = objFolder / f'{obj}_{ext}_stellarFit_Papaderos_{cycle}.png'
        maskPlotFile = objFolder / f'{obj}_{ext}_maskAndFlags_Papaderos_{cycle}.png'
        stellarFluxFile = objFolder / f'{obj}_{ext}_stellarFlux_Papaderos_{cycle}.txt'

        # Load the data
        wave, flux_array, header = sr.import_fits_data(fits_file, instrument='OSIRIS')
        flux = flux_array[idx_band][0] if ext in ('_B', '_R') else flux_array
        lm = sr.LineMesurer(wave, flux, redshift=z_array[i], crop_waves=(wmin_array[i], wmax_array[i]))

        # Starlight wrapper
        sw = SSPsynthesizer()

        # Read output data
        stellar_Wave, obj_input_flux, stellar_flux, fit_output = sw.load_starlight_output(starlight2012_folder/outputFile)
        z_gp = obsData['sample_data']['z_array'][i]
        Mcor, Mint = fit_output['Mcor_tot'], fit_output['Mini_tot']
        mass_galaxy = computeSSP_galaxy_mass(Mcor, 1, z_gp)
        massProcess_galaxy = computeSSP_galaxy_mass(Mint, 1, z_gp)
        idcs_below_20Myr = fit_output['DF'].age_j < 2*10**7
        mass_galaxy_20Myr_percent = np.sum(fit_output['DF'].loc[idcs_below_20Myr, 'Mcor_j'].values)

        # Store starlight configuration values for linux runy
        rc = pn.RedCorr(R_V=RV, E_BV=fit_output['Av_min'] / RV, law=red_law)
        cHbeta_star = rc.cHbetaFromEbv(fit_output['Av_min']/RV)
        starlight_cfg = {'gridFileName': outputFile,
                         'outputFile': outputFile,
                         'saveFolder': starlight2012_folder.as_posix(),
                         'Galaxy_mass_Current': mass_galaxy,
                         'Galaxy_mass_Prosessed': massProcess_galaxy,
                         'Galaxy_mass_Percentbelow20Myr': mass_galaxy_20Myr_percent,
                         'Chi2': fit_output['Chi2'],
                         'A_V_stellarr': fit_output['Av_min'],
                         'cHbeta_stellar': cHbeta_star,
                         'PixelMeanDevPer': fit_output['SumXdev'],
                         'SN': fit_output['SignalToNoise_magnitudeWave']}
        sr.parseConfDict(results_file, starlight_cfg, f'Starlight_run_{cycle}', clear_section=True)

        # Plot the results
        plot_label = f'{obj} spectrum' if ext == '_BR' else f'{obj} blue arm spectrum'
        sw.population_fraction_plots(fit_output, plot_label, 'Mass_fraction', massFracPlotFile, mass_galaxy=mass_galaxy)
        sw.population_fraction_plots(fit_output, plot_label, 'Light_fraction', LightFracPlotFile)
        sw.mask_plot(fit_output, obj, lm.wave, lm.flux, stellar_Wave, stellar_flux, obj_input_flux)#, outputAddress=maskPlotFile)

