from uncertainties import ufloat

a = ufloat(4.702, 0.203)
b = ufloat(2.304, 0.852)
print(a)
print(b)
c = a/b
print(c)
print(c.nominal_value, c.std_dev)