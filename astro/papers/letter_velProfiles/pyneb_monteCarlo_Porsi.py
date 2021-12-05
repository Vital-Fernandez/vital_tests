import numpy as np
import pyneb as pn

def epm_nSII_formula(R_S2, Te):

    a0 = 16.054 - 7.79 / Te - 11.32 * Te
    a1 = -22.66 + 11.08 / Te + 16.02 * Te
    b0 = -21.61 + 11.89 / Te + 14.59 * Te
    b1 = 9.17 - 5.09 / Te - 6.18 * Te
    ne = 1000 * (R_S2 * a0 + a1) / (R_S2 * b0 + b1)

    return ne


print(pn.__version__)
pn.atomicData.includeFitsPath()
pn.atomicData.setDataFileDict("PYNEB_18_01")
print(pn.atomicData.getAllPredefinedDict())

S2 = pn.Atom('S', 2)

# Temperature
TSIII = np.array([12135.0, 188.0])

# Fluxes from Guseva et al 2021
S2_6716A_G = np.array([9.40, 0.31])
S2_6731A_G = np.array([7.77, 0.27])

# Fluxes from Matias
S2_6716A_M = np.array([36.030, 1.916])
S2_6731A_M = np.array([31.056, 2.07])

# Density ratio
RSII_G = S2_6716A_G[0]/S2_6731A_G[0]
RSII_M = S2_6716A_M[0]/S2_6731A_M[0]

# Pyneb densities
neSII_G = S2.getTemDen(RSII_G, tem=TSIII[0], wave1=6716, wave2=6731)
neSII_M = S2.getTemDen(RSII_M, tem=TSIII[0], to_eval='L(6716)/L(6731)')

# Perez-Montero 2017 densities
neSII_G_epm = epm_nSII_formula(RSII_G, TSIII[0])
neSII_M_epm = epm_nSII_formula(RSII_M, TSIII[0])

print(neSII_G_epm)
print()

print('Values from the Guseva et al 2021:')
print(f'- Guseva paper neSII: {225}+/-{76} cm^-3')

print('\nValues calculated with this script')

print(f'\n-Guseva RSII ratio = {RSII_G:.2f}')
print(f'-Guseva Pyneb nSII = {neSII_G:.2f} cm^-3')
print(f'-Guseva EPM-2017 nSII = {neSII_G_epm:.2f} cm^-3')


print(f'\n-Matias RSII ratio = {RSII_M:.2f}')
print(f'-Matias Pyneb nSII = {neSII_M:.2f} cm^-3')
print(f'-Matias EPM-2017 nSII = {neSII_M_epm:.2f} cm^-3')

print(f'\nRatio discrepancy {100 * (1-RSII_G/RSII_M):.1f}%')
print(f'Density discrepancy {100 * (1-neSII_G/neSII_M):.1f}%')
