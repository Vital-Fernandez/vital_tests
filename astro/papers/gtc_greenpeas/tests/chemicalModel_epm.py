import numpy as np
import pyneb as pn
import pandas as pd
import src.specsiser as sr
from astro.papers.gtc_greenpeas.common_methods import compute_arms_flambda

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

columns = ['wavelength', 'flux', 'error', 'int_epm']
indeces = ['O2_3727A', 'O3_4363A', 'H1_4861A', 'O3_4959A', 'O3_5007A', 'O2_7319A', 'S3_6312A', 'S2_6717A', 'S2_6731A']
wavelength = [3727, 4363, 4861, 4959, 5007, 7325, 6312, 6717, 6731]
flux = [5.65E-15, 1.22E-15, 9.14E-15, 1.98E-14, 5.66E-14, 5.49E-16, 1.08E-16, 5.94E-16, 5.45E-16]
eflux = [3.14E-17, 5.61E-18, 9.39E-18, 4.90E-18, 6.33E-18, 2.23E-17, 1.13E-18, 1.43E-17, 1.15E-17]
int_epm = [6276, 1346, 10000, 21657, 61876, 0, 102, 562, 516]


df = pd.DataFrame(index=indeces, columns=columns)
df['wavelength'] = wavelength
df['flux'] = flux
df['eflux'] = eflux
df['int_epm'] = int_epm

cHbeta = 0.017
red_law = obsData['sample_data']['red_law']
RV = obsData['sample_data']['RV']

rc = pn.RedCorr(R_V=RV, law=red_law, cHbeta=cHbeta)
cors = rc.getCorr(df.wavelength.values, rel_wave=4861.0)
lineInts = (df.flux.values / df.loc['H1_4861A'].flux) * cors

f_Hbeta, f_lamdba = compute_arms_flambda(df, red_law, RV, ref_line='H1_4861A')
df['int'] = (df.flux.values / df.loc['H1_4861A'].flux) * np.power(10, cHbeta * f_lamdba)
df['int_epm'] = df['int_epm'] * np.power(10, cHbeta * f_lamdba) / 10000

R_O3 = (df.loc['O3_4959A'].int_epm + df.loc['O3_5007A'].int_epm) / df.loc['O3_4363A'].int_epm
Te_O3_pm2017 = (0.7840 - 0.0001357 * R_O3 + 48.44 / R_O3) * 10000
print('Te_O3_pm2017', Te_O3_pm2017)

Te_S3_pm2017 = (1.19 * (Te_O3_pm2017 / 10000.0) - 0.32) * 10000
print('Te_S3_pm2017', Te_S3_pm2017)

R_S2 = df.loc['S2_6717A'].int_epm / df.loc['S2_6731A'].int_epm
Te = 15618
a0 = 16.054 - 7.79 / Te - 11.32 * Te
a1 = -22.66 + 11.08 / Te + 16.02 * Te
b0 = -21.61 + 11.89 / Te + 14.59 * Te
b1 = 9.17 - 5.09 / Te - 6.18 * Te
ne_S2_pm2017 = 1000 * (R_S2 * a0 + a1) / (R_S2 * b0 + b1)
print('ne_S2_pm2017', ne_S2_pm2017)

R_O2 = df.loc['S2_6717A'].int_epm / df.loc['S2_6731A'].int_epm
ne = 412.0
Te = 15618.0
TO2_Hag2006 = ((1.2 + 0.002 * ne + 4.2/ne) / ((1/(Te/10000.0)) + 0.08 + 0.003*ne + 2.5/ne))*10000.0
print('TO2_Hag2006', TO2_Hag2006)

O3_coeff = (df.loc['O3_4959A', 'int_epm'] + df.loc['O3_5007A', 'int_epm']) / df.loc['H1_4861A', 'int_epm']
O3_abund = np.log10(O3_coeff) + 6.1868 + 1.2491/(Te/10000.0) - 0.5816 * np.log10(Te/10000.0)
print('O3_abund', O3_abund)

T_low = TO2_Hag2006
S2_coeff = (df.loc['S2_6717A', 'int_epm'] + df.loc['S2_6731A', 'int_epm']) / df.loc['H1_4861A', 'int_epm']
S2_abund = np.log10(S2_coeff) + 5.463 + 0.941/(T_low/10000.0) - 0.97 * np.log10(T_low/10000.0)
print('S2_abund', S2_abund)
