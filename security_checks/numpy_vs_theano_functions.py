import theano
import numpy as np
import src.specsiser as sr
import theano.tensor as T
import pyneb as pn

theano.config.on_unused_input = 'ignore'

print('\nTheano example')
x, y = T.dscalar('x'), T.dscalar('y')
z = x + y
print('z = ', z.eval({x : 16.3, y : 12.1}))

print('\nFlux computation using numpy arrays')
label_list = ['O3_5007A', 'S3_6312A', 'H1_6563A', 'He1_5876A']
ion_list = ['O3', 'S3', 'H1r', 'He1r']

emtt = sr.EmissionTensors(label_list, ion_list)

emis_ratio, cHbeta, flambda, abund, ftau = 0.352, 0.12, 0.2, 7.8, 0.0
params = dict(emis_ratio=emis_ratio, cHbeta=cHbeta, flambda=flambda, abund=abund, ftau=ftau)
flux_i = emtt.emFluxEqDict['O3_5007A'](emis_ratio, cHbeta, flambda, abund, ftau, {})
print('Flux_i', flux_i)

print('\nFlux computation using theano graphs')
emis_ratio_t, cHbeta_t, flambda_t, abund_t, ftau_t = T.dscalars('emis_ratio_t', 'cHbeta_t', 'flambda_t', 'abund_t', 'ftau_t')
emGraph = sr.EmissionFluxModel(label_list, ion_list)
params = {emis_ratio_t: 0.352, cHbeta_t: 0.12, flambda_t: 0.2, abund_t: 7.8, ftau_t: 0.0}
flux_j = emGraph.emFluxEqDict['O3_5007A'](emis_ratio_t, cHbeta_t, flambda_t, abund_t, ftau_t, {})
print('flux_j = ', flux_j.eval(params))


print('\nEmissivity interpolation')



