import pyneb as pn
import numpy as np
import src.specsiser as sr


def corO2_7319(emis_ratio, cHbeta, flambda, O2_abund, O3_abund, Te_high):
    fluxCorr = np.power(10, O2_abund + emis_ratio - flambda * cHbeta - 12) + \
               np.power(10, O3_abund + 0.9712758487381 + np.log10(np.power(Te_high/10000.0, 0.44)) - flambda * cHbeta - 12)

    return fluxCorr


Te, ne = 10000.0, 100.0
cHbeta = 0.255
H1_ion = pn.RecAtom('H', 1)
O2_ion = pn.Atom('O', 2)
O3_ion = pn.Atom('O', 3)
cHbeta = 0.12
flambda = - 0.55

O2_abund = 0.00000648
O3_abund = 0.00004578

logO2abund = 12 + np.log10(O2_abund)
logO3abund = 12 + np.log10(O3_abund)

Hbeta_emis = H1_ion.getEmissivity(tem=Te, den=ne, wave=4861)
O2_emis = (O2_ion.getEmissivity(tem=Te, den=ne, wave=7319) + O2_ion.getEmissivity(tem=Te, den=ne, wave=7319))/Hbeta_emis
logO2_emis = np.log10(O2_emis)

f_O2 = O2_abund * O2_emis * np.power(10, -cHbeta * flambda)
flog_O2 = np.power(10, logO2abund + logO2_emis - flambda * cHbeta - 12)

f_O2_rec = (O3_abund * 9.36 * np.power(Te/10000.0, 0.44)) * np.power(10, -cHbeta * flambda)
flog_O2_rec = np.power(10, logO3abund + np.log10(9.36) + np.log10(np.power(Te/10000.0, 0.44)) - flambda * cHbeta - 12)

f_O2_tot = f_O2 + f_O2_rec
flog_O2_tot = flog_O2 + flog_O2_rec

emtt = sr.EmissionTensors(['O2_7319A_b', 'O3_5007A'], ['O2', 'O3'])
f_O2_tt = emtt.compute_flux('O2_7319A_b',
                             logO2_emis,
                             cHbeta,
                             flambda,
                             logO2abund,
                             0.0,
                             O3=logO3abund,
                             T_high=Te)


# print(f_O2)
# print(flog_O2)
# print
# print(f_O2_rec)
# print(flog_O2_rec)
# print
print(f_O2_tot)
print(flog_O2_tot)
print(np.log10(f_O2_tot))
print(np.log10(flog_O2_tot))
print(f_O2_tt)

# print(emisTT.emFluxTensors['O2_7319A_b'](np.log10(O2_emis), cHbeta, flambda, logO2abund, logO3abund, temp))
# print(corO2_7319(np.log10(O2_emis), cHbeta, flambda, logO2abund, logO3abund, temp)

