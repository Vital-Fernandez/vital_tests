import pyneb as pn
import numpy as np

redCor1 = pn.RedCorr(E_BV=0.193, R_V=3.1, law='CCM89')
redCor2 = pn.RedCorr(E_BV=0.134, R_V=3.1, law='CCM89')
redCor3 = pn.RedCorr(E_BV=0.170, R_V=3.1, law='CCM89')
redCor4 = pn.RedCorr(E_BV=0.350, R_V=3.1, law='CCM89')

ext_array = np.array([redCor1.cHbeta, redCor2.cHbeta, redCor3.cHbeta])
den_array = np.array([140, 120, 100])
O_array = np.array([7.81, 7.82, 7.74])

print(redCor1.cHbeta, redCor2.cHbeta, redCor3.cHbeta, redCor4.cHbeta)
print(np.mean(ext_array))
print(np.mean(den_array))
print(np.mean(O_array))

def undo(twuelf_XH):
    return np.power(10, twuelf_XH-12)

O_ratio = undo(7.71)/undo(7.36)
S_ratio = undo(5.91)/undo(5.26)
eta = O_ratio/S_ratio

print('-------')
print(O_ratio)
print(S_ratio)
print(eta)
print(np.log10(eta))
