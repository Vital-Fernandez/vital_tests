import numpy as np
import pandas as pd
from pathlib import Path
import src.specsiser as sr
from src.specsiser.treatment import fits_db
from astro.data.muse.common_methods import grid_columns
from timeit import default_timer as timer

# Declare data and files location
obsData = sr.loadConfData('../muse_greenpeas.ini', group_variables=False)
objList = obsData['data_location']['object_list']
fileList = obsData['data_location']['file_list']
fitsFolder = Path(obsData['data_location']['fits_folder'])
dataFolder = Path(obsData['data_location']['data_folder'])
resultsFolder = Path(obsData['data_location']['results_folder'])

z_objs = obsData['sample_data']['z_array']
pertil_array = obsData['sample_data']['percentil_array']
noise_region = obsData['sample_data']['noiseRegion_array']
norm_flux = obsData['sample_data']['norm_flux']

# Load the data
file_address = Path('/home/vital/Dropbox/Astrophysics/Data/muse-Ricardo/Data/HII-CHI-mistry_1Myr_grid_O.txt')
grid_3D_DF = pd.read_csv(file_address, delim_whitespace=True, header=0)

# Prepare interpolators
model_variables = ['logOH', 'logU', 'logNO']
gw = sr.ModelGridWrapper()
grid_dict, axes_cords_a = gw.ndarray_from_DF(grid_3D_DF, axes_columns=model_variables)
grid_interpolators = gw.generate_xo_interpolators(grid_dict, model_variables, axes_cords_a, interp_type='point')

model_lines = np.array(list(grid_dict.keys()))

lines_MUSE = np.array(['O3_4959A', 'O3_5007A',
                       'He1_5876A',
                       'N2_6548A', 'H1_6563A', 'N2_6584A',
                       'S2_6716A', 'S2_6731A',
                       'S3_9069A', 'S3_9531A', 'R_O3N2', 'R_S2S3'])

ratio_format = {'R_S2S3': ['-', '-', '$R_{S2S3}$'],
                'R_O3N2': ['-', '-', '$R_{O3N2}$']}

idcs_lines = np.isin(lines_MUSE, model_lines)
input_lines = lines_MUSE[idcs_lines]

logOH_range = np.round(np.linspace(7.15, 9.0, 5), 3)
logU_range = np.round(np.linspace(-3.90, -1.40, 5), 3)
logNO_range = np.round(np.linspace(-1.90, -0.01, 5), 3)

# Default configuration
conf_file = Path(dataFolder/'CGCG007/chemical_fitting_conf.txt')
conf_fit = sr.loadConfData(conf_file)

objFolder = resultsFolder / f'CGCG007/'
outputFits = objFolder / f'grid_logOH_logU_logNO_advi.fits'

# Loop throught the grid of synthetic conditions (around 30 seconds per fit)
start = timer()
i_step, n_steps = 0, logNO_range.size * logOH_range.size * logNO_range.size
for i, logOH in enumerate(logOH_range):
    for j, logU in enumerate(logU_range):
        for k, logNO in enumerate(logNO_range):

            i_step += 1
            if i_step > 59:
                print(f'- Fit {i_step}/{n_steps}: expected time {n_steps*30/3600:0.3f} hours')

                # True value coordinate for interpolation
                coord_true = [[logOH, logU, logNO]]
                header_params = {'logOH': logOH, 'logU': logU, 'logNO': logNO}

                # Output files
                cord_label = f'{logOH*1000:.0f}_{logU*-1000:.0f}_{logNO*-1000:.0f}'
                outputCfg = objFolder/f'{cord_label}.txt'
                outputDb = objFolder/f'{cord_label}'

                # Fill the dataframe with integrated flux
                linesDF = pd.DataFrame(index=input_lines, columns=['wavelength', 'ion',	'intg_flux', 'intg_err'])
                for line in input_lines:
                    if not line.startswith('R_'):
                        ion, wavelength, latexLabel = sr.label_decomposition(line, scalar_output=True)
                        flux = np.power(10, grid_interpolators[line](coord_true).eval()[0][0])
                        linesDF.loc[line, :] = wavelength, ion, flux, flux * 0.05
                    else:
                        ion, wavelength, latexLabel = 'None', 0.0, line
                        flux = np.power(10, grid_interpolators[line](coord_true).eval()[0][0])
                        linesDF.loc[line, :] = wavelength, ion, flux, flux * 0.05

                # Declare sampler
                obj1_model = sr.SpectraSynthesizer()

                # Declare region physical model
                obj1_model.define_region(linesDF)

                # Declare region physical model
                obj1_model.simulation_configuration(model_parameters=conf_fit['inference_model_configuration']['parameter_list'],
                                                    prior_conf_dict=conf_fit['priors_configuration'],
                                                    grid_interpolator=grid_interpolators)

                obj1_model.photoionization_sampling(conf_fit['inference_model_configuration']['parameter_list'])

                obj1_model.run_sampler(1000, 3000, nchains=8, njobs=4)

                obj1_model.save_fit(outputFits, cord_label, output_format='fits', user_header=header_params)

end = timer()
print(f'Working time:{(end-start)/n_steps} seconds per fit')
