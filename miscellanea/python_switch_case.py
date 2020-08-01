# def opA_2(x, y):
#     return x - y
#
#
# def opB_3(x, y, z):
#     return x + y - z
#
#
# def opC_2(x, y):
#     return x * y
#
#
# def opD_3(x, y, z):
#     return x * y + z
#
#
# def dispact_func(op_dict, op, x, y, **kwargs):
#     return op_dict.get(op, lambda x, y, kwargs: None)(x, y, kwargs)
#
#
# operation_dict = {
#                 'opA_2': opA_2,
#                 'opB_3': opB_3,
#                 'opC_2': opC_2,
#                 'opD_3': opD_3
#                 }
#
# op_dict_lambda = {}
# for label, func in operation_dict.items():
#     if '2' in label:
#         op_dict_lambda[label] = lambda x, y, kwargs: func(x, y)
#     if '3' in label:
#         op_dict_lambda[label] = lambda x, y, kwargs: func(x, y, kwargs['z'])
#
# coefs_dict = {'i': 1, 'j': 2, 'k': 3, 'z': 4}
#
# print('True result ', operation_dict['opB_3'](1, 2, coefs_dict['z']))
# print('dispact_func result ', dispact_func(op_dict_lambda, 'opB_3', 1, 2, **coefs_dict))



def opA_2(x, y):
    return x + y


def opB_3(x, y, z):
    return x + y - z


def opC_2(x, y):
    return x * y


def opD_3(x, y, z):
    return x * y + z


op_dict = {'opA_2': opA_2,
           'opB_3': opB_3,
           'opC_2': opC_2,
           'opD_3': opD_3
           }

op_lambda_dict = {'opA_2': lambda x, y, kwargs: op_dict['opA_2'](x, y),
                  'opB_3': lambda x, y, kwargs: op_dict['opB_3'](x, y, kwargs['z']),
                  'opC_2': lambda x, y, kwargs: op_dict['opC_2'](x, y),
                  'opD_3': lambda x, y, kwargs: op_dict['opD_3'](x, y, kwargs['z']),
                  }


def dispatch_op(func_dict, op, x, y, **kwargs):
    return func_dict.get(op, lambda a, b, c: None)(x, y, kwargs)

coefs_dict = {'i': 1, 'j': 2, 'k': 3, 'z': 4}

print('Original lambda dict result:', dispatch_op(op_lambda_dict, 'opD_3', 1, 2, **coefs_dict))

# op_looplambda_dict = {}
#
#
# def add_to_dict(function, my_dict, my_key):
#     if '2' in label:
#         my_dict[my_key] = lambda x, y, kwargs: function(x, y)
#     if '3' in label:
#         my_dict[my_key] = lambda x, y, kwargs: function(x, y, kwargs['z'])
#
#
# for label, func in op_dict.items():
#     add_to_dict(func, op_looplambda_dict, label)
#
#
# print('Loop lambda dict result:', dispatch_op(op_looplambda_dict, 'opD_3', 1, 2, **coefs_dict))

op_looplambda_dict = {}


def add_to_dict(function, my_dict, my_key):
    if '2' in label:
        my_dict[my_key] = lambda x, y, kwargs: function(x, y)
    if '3' in label:
        my_dict[my_key] = lambda x, y, kwargs: function(x, y, kwargs['z'])


for label, func in op_dict.items():
    add_to_dict(func, op_looplambda_dict, label)


print('Loop lambda dict result:', dispatch_op(op_looplambda_dict, 'opD_3', 1, 2, **coefs_dict))


# op_dict_lambda = {}
# for label, func in operation_dict.items():
#     if '2' in label:
#         op_dict_lambda[label] = lambda x, y, kwargs: func(x, y)
#     if '3' in label:
#         op_dict_lambda[label] = lambda x, y, kwargs: func(x, y, kwargs['z'])



# import theano.tensor as tt
# from theano import function
#
# def sum2(x, y):
#     return x + y
#
#
# def sum3(x, y, z):
#     return x + y + z
#
#
# def mul2(x, y):
#     return x * y
#
#
# def mul3(x, y, z):
#     return x * y * z
#
#
#
# op_dict = {
#     'sum2': sum2,
#     'sum3': sum3,
#     'mul2': mul2,
#     'mul3': mul3
#     }
#
# paramDict = {
#             'sum2': ['x', 'y'],
#             'sum3': ['x', 'y', 'z'],
#             'mul2': ['x', 'y'],
#             'mul3': ['x', 'y', 'z']
#              }
#
#
# op_dict_new = {}
#
# for label, func in op_dict.items():
#     func_params = tt.dscalars(paramDict[label])
#     if '2' in label:
#         op_dict_new[label] = lambda x=0, y=0, kwargs={}: func(x,y)
#     if '3' in label:
#         op_dict_new[label] = lambda x=0, y=0, kwargs={}: func(x, y, kwargs['z'])
#
#
# op_dict_tt = {}
#
# for label, func in op_dict.items():
#     func_params = tt.dscalars(paramDict[label])
#     if '2' in label:
#         op_dict_tt[label] = lambda x=0, y=0, kwargs={}: func(x, y)
#     if '3' in label:
#         op_dict_tt[label] = lambda x=0, y=0, kwargs={}: func(x, y, kwargs['z'])
#
#
# def dispatch_dict(op_dict, op, x, y, **kwargs):
#     return op_dict.get(op, lambda x, y, kwargs: None)(x, y, kwargs)
#
#
# myDict = {'i': 1, 'j': 2, 'k': 3, 'z': 4}
# # print(dispatch_dict(op_dict, 'mul3', 1, 2, **myDict))
# print(dispatch_dict(op_dict_new, 'sum3', 1, 2, **myDict))
# print(dispatch_dict(op_dict_tt, 'sum3', 1, 2, **myDict))


# myDict = {'a':1, 'b':2, 'c':3, 'd': 6}
# myDict2 = {'a':1, 'b':2, 'c':3, 'd': 7}
#
# print(myDict)
# for key, value in myDict.items():
#     myDict[key] = value * 2
# print(myDict)
#

# print(myDict)
# print(myDict2)
# myDict2['d'] = myDict['d']
# print(myDict)
# print(myDict2)
# myDict['d'] = 10
# print(myDict)
# print(myDict2)


# import numpy as np
# from src.specsiser.physical_model.gasEmission_functions import EmissionFluxModel, EmissionTensors
# from timeit import default_timer as timer
#
# lineLabels = np.array(['H1_4341A', 'O3_4363A', 'He1_4471A', 'He2_4686A', 'Ar4_4740A',
#                        'O3_4959A', 'O3_5007A', 'He1_5876A', 'S3_6312A', 'N2_6548A',
#                        'H1_6563A', 'N2_6584A', 'He1_6678A', 'S2_6716A', 'S2_6731A',
#                        'He1_7065A', 'Ar3_7136A', 'O2_7319A', 'O2_7319A_b', 'O2_7330A',
#                        'S3_9069A', 'S3_9531A'])
#
# lineIons = np.array(['H1r', 'O3', 'He1r', 'He2r', 'Ar4', 'O3', 'O3', 'He1r', 'S3', 'N2',
#                      'H1r', 'N2', 'He1r', 'S2', 'S2', 'He1r', 'Ar3', 'O2', 'O2', 'O2',
#                      'S3', 'S3'])
#
# emEq = EmissionTensors(lineLabels, lineIons)
#
# lineLabel, lineIon = 'O3_5007A', 'O3'
#
# flux = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.0)
#
# print(flux)
#
# start = timer()
# fluxa = emEq.compute_flux(lineLabel, 2.0, 0.60, 0.12, 7, 0.14)
# end = timer()
# totalTime = end - start
# print('Time 1 000 000 (from 1)', totalTime * 1000000)
#
# start = timer()
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# fluxb = emEq.compute_flux(lineLabel, 2.0, 0.58, 0.2, 7, 0.0)
# fluxc = emEq.compute_flux(lineLabel, 3.0, 0.08, 0.2, 6, 0.2)
# fluxd = emEq.compute_flux(lineLabel, 1.0, 0.58, 0.2, 5, 0.4)
# fluxe = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 4, 0.2)
# fluxa = emEq.compute_flux(lineLabel, 3.0, 0.58, 0.2, 8, 0.1)
# end = timer()
# totalTime = end - start
# print('Time 1 000 000 (from 100)', (totalTime/100.0) * 1000000)
#
#
# def test_theano_function_kwargs():
#     import numpy as np
#     f = theano_function([x, y, z], [x+y], dim=1, on_unused_input='ignore')
#     assert np.linalg.norm(f([1, 2], [3, 4], [0, 0]) - np.asarray([4, 6])) < 1e-9
#
#     f = theano_function([x, y, z], [x+y], dtypes={x: 'float64', y: 'float64'},
#                                      dim=1, on_unused_input='ignore')
#     xx = np.arange(3).astype('float64')
#     yy = 2*np.arange(3).astype('float64')
#     zz = 2*np.arange(3).astype('float64')
#     assert np.linalg.norm(f(xx, yy, zz) - 3*np.arange(3)) < 1e-9





# def sum2(x, y):
#     return x + y
#
#
# def sum3(x, y, z):
#     return x + y + z
#
#
# def mul2(x, y):
#     return x * y
#
#
# def mul3(x, y, z):
#     return x * y * z
#
#
# myparams = ['x', 'y', 'kwargs']
# op_dict = {
#     'sum2': lambda x, y, kwargs: sum2(x, y),
#     'sum3': lambda x, y, kwargs: sum3(x, y, kwargs['z']),
#     'mul2': lambda x, y, kwargs: mul2(x, y),
#     'mul3': lambda x, y, kwargs: mul3(x, y, kwargs['z']),
# }
#
#
# def dispatch_dict(op_dict, op, x, y, **kwargs):
#     return op_dict.get(op, lambda a, b, c: None)(x, y, kwargs)
#
#
# myDict = {'i': 1, 'j': 2, 'k': 3, 'z': 4}
#
# a = dispatch_dict(op_dict, 'sum3', 1, 2, **myDict)
# print(a)
#


# import numpy as np
# from src.specsiser.physical_model.gasEmission_functions import EmissionTensors
#
# emtt = EmissionTensors()

# O2_abund = 0.5*1e-7
# O2_emiss = 0.582
# cHbeta = 0.355
#
# O3_abund = 0.5*1e-8
# O3_emiss = 0.758
#
# flambda = 0.150
#
# O2_flux = emtt.fluxEq['metals'](O2_emiss, cHbeta, flambda, O2_abund, 0.0)
# O3_flux = emtt.fluxEq['metals'](O3_emiss, cHbeta, flambda, O3_abund, 0.0)
#
# O2_flux_log = emtt.fluxEq_log['metals'](np.log10(O2_emiss), cHbeta, flambda, 12 + np.log10(O2_abund), 0.0)
# O3_flux_log = emtt.fluxEq_log['metals'](np.log10(O3_emiss), cHbeta, flambda, 12 + np.log10(O3_abund), 0.0)
#
# O2_O3_log = emtt.fluxEq_log['O2_7319A_b'](np.log10(O2_emiss), cHbeta, flambda, 12 + np.log10(O2_abund), 0.0,
#                                           12 + np.log10(O3_abund), np.log10(O3_emiss))
#
# print('Flux O2', O2_flux)
# print('Flux O3', O3_flux)
# print('Flux O2 log', np.power(10, O2_flux_log))
# print('Flux O3 log', np.power(10, O3_flux_log))
# print('Flux sum', O2_flux + O3_flux)
# print('Flux sum log', np.power(10, O2_flux_log) + np.power(10, O3_flux_log))
# print('Flux sumlog', np.log10(O2_flux + O3_flux))
# print('Flux O2 log', O2_flux_log)
# print('Flux O3 log', O3_flux_log)
# print('Flux log(O2+O3)', np.log10(O2_flux + O3_flux))
# print('Flux sum', np.log10(np.power(10, O2_flux_log) + np.power(10, O3_flux_log)))
# print('O2_O3_log', O2_O3_log)


# def sum2(x, y):
#     return x + y
#
#
# def sum3(x, y, z):
#     return x + y + z
#
#
# def mul2(x, y):
#     return x * y
#
#
# def mul3(x, y, z):
#     return x * y * z
#
#
# def dispatch_dict(op, x, y, **kwargs):
#     return {
#         'sum2': lambda: sum2(x, y),
#         'sum3': lambda: sum3(x, y, kwargs['z']),
#         'mul2': lambda: mul2(x, y),
#         'mul3': lambda: mul3(x, y, kwargs['z'])
#     }.get(op, lambda: None)()
#
#
# print(dispatch_dict('sum2', 1, 2, z=3))
# print(dispatch_dict('sum3', 1, 2, z=3))
#
# print(dispatch_dict('sum3', 1, 2, z=3))


# def sum(x, y):
#     return x + y
#
#
# def sub(x, y):
#     return x - y
#
#
# def mul(x, y):
#     return x * y
#
#
# def div(x, y):
#     return x / y
#
#
# operDict = {
#         'add': sum,
#         'sub': sub,
#         'mul': mul,
#         'div': div}
#
# operlambdaDict = {
#         'add': lambda x, y: x + y,
#         'sub': lambda x, y: x - y,
#         'mul': lambda x, y: x * y,
#         'div': lambda x, y, z: x / y,
# }
#
# operlambdaFuncDict = {
#     'sum': lambda: sum,
#     'sub': lambda: sub,
#     'mul': lambda: mul,
#     'div': lambda: div
#     }
#
#
#
# def dispatch_dict7(op, x, y):
#     return {
#         'sum': lambda: sum(x, y),
#         'sub': lambda: sub(x, y),
#         'mul': lambda: mul(x, y),
#         'div': lambda: div(x, y)
#     }.get(op, lambda: None)()
#
#
# def dispatch_dict(operator, x, y):
#     return {
#         'add': lambda: x + y,
#         'sub': lambda: x - y,
#         'mul': lambda: x * y,
#         'div': lambda: x / y,
#     }.get(operator, lambda: None)()
#
#
# def dispatch_dict2(operator, x, y):
#     return {
#         'add': lambda: sum(x, y),
#         'sub': lambda: sub(x, y),
#         'mul': lambda: mul(x, y),
#         'div': lambda: div(x, y),
#     }.get(operator, lambda: None)()
#
#
# def dispatch_dict3(operator, x, y):
#     return operDict.get(operator, lambda: None)(x, y)
#
#
# def dispatch_dict4(operator, **kwargs):
#     return operDict.get(operator, None)(**kwargs)
#
#
# def dispatch_dict5(operator, x, **kwargs):
#     return operDict.get(operator, None)(x, **kwargs)
#
#
# def dispatch_dict6(operator, x, **kwargs):
#     return operlambdaDict.get(operator, None)(x, **kwargs)
#
#
#
#
#
#
#
# myDict = {'y': 10, 'z': 20}
#
# suml = lambda x, y: x+y
# print(suml(0, 10))
# print(dispatch_dict('add', 10, 10))
# print(dispatch_dict2('add', 20, 10))
# print(dispatch_dict3('add', 30, 10))
# print(dispatch_dict4('add', x=40, y=10))
# print(dispatch_dict5('add', 50, y=10))
# print(dispatch_dict6('add', 60, y=10))
# print(dispatch_dict6('add', 70, **myDict))

# myDict = {'abund': 10.0, 'ftau': 0.0, 'emis_ratio': 10.0}
#
# print('H1r', emtt.fluxEq_log['H1r'](30, 10, 2, None, None))
# print('Metals', emtt.fluxtt['metals'](30, 10, 2, 10, 0.0))
# print('O2_7319A_b', emtt.fluxtt['O2_7319A_b'](30, 10, 2, 10, 0.0, 10, 2))
# print('Metals 2', emtt.fluxtt['metals'](30, 10, 2, **myDict))

# print('Coso', emtt.fluxtt['O2_7319A_b'](30, 10, 2, None, None, 30, 10))
