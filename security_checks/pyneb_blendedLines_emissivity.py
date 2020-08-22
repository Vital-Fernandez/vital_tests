import numpy as np
import pyneb as pn

H1 = pn.RecAtom('H', 1)
O2 = pn.Atom('O', 2)
S2 = pn.Atom('S', 2)

S2.printIonic()


ne, Te, abund = 150.0, 12000.0, 0.0005

Hbeta_em = H1.getEmissivity(Te, ne, wave=4861)
O2_3726A_em = O2.getEmissivity(Te, ne, wave=3726)
O2_3729A_em = O2.getEmissivity(Te, ne, wave=3729)
O2_3720b_em = O2_3726A_em + O2_3729A_em

O2_3726A_flux = abund * O2_3726A_em/Hbeta_em
O2_3729A_flux = abund * O2_3729A_em/Hbeta_em
O2_3720b_flux = abund * O2_3720b_em/Hbeta_em

O2_3726A_ref = 'I(3,1)'
O2_3729A_ref = 'I(2,1)'
O2_3720b_ref = 'I(2,1)+I(3,1)'

O2_3726A_abund = O2.getIonAbundance(O2_3726A_flux, Te, ne, to_eval=O2_3726A_ref, Hbeta=1)
O2_3729A_abund = O2.getIonAbundance(O2_3729A_flux, Te, ne, to_eval=O2_3729A_ref, Hbeta=1)
O2_3720b_abund = O2.getIonAbundance(O2_3720b_flux, Te, ne, to_eval=O2_3720b_ref, Hbeta=1)

print('O2_3726A_abund', O2_3726A_abund)
print('O2_3729A_abund', O2_3729A_abund)
print('O2_3720b_abund', O2_3729A_abund)
