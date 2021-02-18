import scipy as spy
import exoplanet as xo
import numpy as np
import pymc3
import theano.tensor as tt
from src.specsiser.inference_model import displaySimulationData
from astro.papers.gtc_greenpeas.common_methods import reading_epm_grids
import arviz as az
import matplotlib.pyplot as plt

grid_line_dict, params_range_dict = reading_epm_grids()

print(params_range_dict)

grid_line_dict.pop('S3_9069A')

# Generate interpolators
lineInterpolators = {}
for lineLabel, lineGrid in grid_line_dict.items():
    xo_interp = xo.interp.RegularGridInterpolator([params_range_dict['logU'],
                                                 params_range_dict['Teff'],
                                                 params_range_dict['OH']], lineGrid)
    lineInterpolators[lineLabel] = xo_interp.evaluate

# # ------------------------------------------- Testing
# lineLabel_t = 'O2_3726A_m'
# temp_t, logU_t, OH_t = 55000, -2.0, 7.4
# lineCube = grid_line_dict[lineLabel_t]
#
# # Grid point
# idx_temp = params_range_dict['Teff'] == temp_t
# idx_logU = params_range_dict['logU'] == logU_t
# idx_OH = params_range_dict['OH'] == OH_t
# print('True value', lineCube[idx_logU, idx_temp, idx_OH])
#
# # Scipy interpolation
# spy_RGridInterp = spy.interpolate.RegularGridInterpolator((params_range_dict['logU'],
#                                                            params_range_dict['Teff'],
#                                                            params_range_dict['OH']), lineCube)
# print('Scipy interpolation', spy_RGridInterp([[logU_t, temp_t, OH_t]]))
#
# # Exoplanet interpolation
# exop_interp = lineInterpolators[lineLabel_t]
# coordB = np.stack(([logU_t], [temp_t], [OH_t]), axis=-1)
# print('Exoplanet interpolation', exop_interp(coordB).eval())
#
# Creat synthetic observations
temp_true, logU_true, OH_true = 53500.0, -1.70, 7.4
temp_true, logU_true, OH_true = 48620.0, -2.15, 7.4
coord_true = np.stack(([logU_true], [temp_true], [OH_true]), axis=-1)
lineFluxes = np.zeros(len(grid_line_dict.keys()))
for i, item in enumerate(lineInterpolators.items()):
    lineLabel, lineInterpolator = item
    lineFluxes[i] = lineInterpolator(coord_true).eval()[0][0]
lineErr = lineFluxes * 0.05
print(lineFluxes)
print(lineErr)

# Inference model
lineLabels = np.array(list(grid_line_dict.keys()))
lineRange = np.arange(lineLabels.size)

with pymc3.Model() as model:

    # Priors
    OH = OH_true
    Teff = pymc3.Uniform('Teff', lower=30000.0, upper=90000.0)
    logU = pymc3.Uniform('logU', lower=-4, upper=-1.5)

    # Interpolation coord
    grid_coord = tt.stack([[logU], [Teff], [OH]], axis=-1)

    # Loop throught
    for i in lineRange:

        # Line intensity
        lineInt = lineInterpolators[lineLabels[i]](grid_coord)

        # Inference
        Y_emision = pymc3.Normal(lineLabels[i], mu=lineInt, sd=lineErr[i], observed=lineFluxes[i])

    displaySimulationData(model)

    trace = pymc3.sample(5000, tune=2000, chains=2, cores=1, model=model)

print(trace)
print(pymc3.summary(trace))

az.plot_trace(trace)
plt.show()
az.plot_posterior(trace)
plt.show()

# Generate interpolators