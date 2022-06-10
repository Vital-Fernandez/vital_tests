import numpy as np
from uncertainties import ufloat
from uncertainties import umath
# a = ufloat(4.702, 0.203)
# b = ufloat(2.304, 0.852)
# print(a)
# print(b)
# c = a/b
# print(c)
# print(c.nominal_value, c.std_dev)

a_value = 10.50
b_value = 2.00
a_error = a_value * 0.02
b_error = b_value * 0.02

a = ufloat(a_value, a_value * 0.02)
b = ufloat(b_value, b_value * 0.02)

print(a)
print(b)

a_log = umath.log10(a)
b_log = umath.log10(b)

print(a_log)
print(b_log)
print('Umath: a ->', a_log.nominal_value, a_log.std_dev)
print('Umath: b ->', b_log.nominal_value, b_log.std_dev)
print(np.log10(a))
print(np.log10(b))

# a_conv = umath.pow(10, a_log)
# b_conv = umath.pow(10, b_log)
#
# print(a_conv)
# print(b_conv)
#
# a_mine = [np.log10(a_value), np.log10(1.0 + (a_value*0.02)/a_value)]
# b_mine = [np.log10(b_value), np.log10(1.0 + (b_value*0.02)/b_value)]
#
# print('a_mine', a_mine)
# print('b_mine', b_mine)
#
# a_form1 = [np.log10(a_value), 0.434 * (a_value*0.02)/a_value]
# b_form1 = [np.log10(b_value), 0.434 * (b_value*0.02)/b_value]
#
# print(a_form1)
# print(b_form1)
#
# a_form2 = [np.log10(a_value), (a_value*0.02)/(a_value*np.log(10))]
# b_form2 = [np.log10(b_value), (b_value*0.02)/(b_value*np.log(10))]
#
# print(a_form2)
# print(b_form2)
#
# log10_factor = 1.0 / np.log(10)
# print('Factor log10: ', 1.0/np.log(10))
#
# a_conv = [np.log10(a_value), a_error/a_value * log10_factor]
# b_conv = [np.log10(b_value), b_error/b_value * log10_factor]
#
# print(a_conv)
# print(b_conv)

# inputGridFlux = np.log10(self.emissionFluxes.astype(float))
# inputGridErr = np.log10(1.0 + self.emissionErr.astype(float) / self.emissionFluxes.astype(float))



