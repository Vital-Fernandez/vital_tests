import numpy as np
import pyneb as pn
import pandas as pd
import src.specsiser as sr
from astro.papers.gtc_greenpeas.common_methods import compute_arms_flambda
from src.specsiser.components.chemical_model import Standard_DirectMetchod
import pyneb as pn

conf_file_address = '../../../papers/gtc_greenpeas/gtc_greenpeas_data.ini'
obsData = sr.loadConfData(conf_file_address, objList_check=True, group_variables=False)

obj = 'gp121903'
columns = ['wavelength', 'flux', 'eflux', 'int_epm']
indeces = ['O2_3727A', 'O3_4363A', 'H1_4861A', 'O3_4959A', 'O3_5007A', 'N2_6584A', 'S3_6312A', 'S2_6716A', 'S2_6731A']
wavelength = np.array([3727, 4363, 4861, 4959, 5007, 6584, 6312, 6717, 6731])
flux = np.array([5.65E-15, 1.22E-15, 9.14E-15, 1.98E-14, 5.66E-14, 1.366E-015, 1.08E-16, 5.94E-16, 5.45E-16])
eflux = np.array([3.14E-17, 5.61E-18, 9.39E-18, 4.90E-18, 6.33E-18, 1.11E-017, 1.13E-18, 1.43E-17, 1.15E-17])
int_epm = np.array([6276, 1346, 10000, 21657, 61876, 1297, 102, 562, 516])
int_err = np.array([40, 7, 20, 26, 69, 12, 1, 23, 13])

df = pd.DataFrame(index=indeces, columns=columns)
df['wavelength'] = wavelength
df['flux'] = flux
df['eflux'] = eflux
df['int_epm'] = int_epm
df['int_err'] = int_err
df['rel_flux'] = df.flux.values / df.loc['H1_4861A'].flux * 10000

R_O3 = (df.loc['O3_4959A'].int_epm + df.loc['O3_5007A'].int_epm) / df.loc['O3_4363A'].int_epm
Te_O3_pm2017 = (0.7840 - 0.0001357 * R_O3 + 48.44 / R_O3) * 10000
print('\nTe_O3_pm2017', Te_O3_pm2017)

Te = Te_O3_pm2017/10000
R_S2 = df.loc['S2_6731A'].int_epm / df.loc['S2_6731A'].int_epm
a0 = 16.054 - 7.79 / Te - 11.32 * Te
a1 = -22.66 + 11.08 / Te + 16.02 * Te
b0 = -21.61 + 11.89 / Te + 14.59 * Te
b1 = 9.17 - 5.09 / Te - 6.18 * Te
ne_S2_pm2017 = 1000 * (R_S2 * a0 + a1) / (R_S2 * b0 + b1)
print('\nne_S2_pm2017', ne_S2_pm2017)

# Pyneb Electron parameters
O2, O3 = pn.Atom('O', 2), pn.Atom('O', 3)
S2, S3 = pn.Atom('S', 2), pn.Atom('S', 3)
N2 = pn.Atom('N', 2)

neSII = S2.getTemDen(df.loc['S2_6716A', 'int_epm']/df.loc['S2_6731A', 'int_epm'], tem=15618, wave1=6716, wave2=6731)
TeOIII = O3.getTemDen(df.loc['O3_4363A', 'int_epm']/df.loc['O3_5007A', 'int_epm'], den=412, wave1=4363, wave2=5007)
print(f'\nPyneb getTemDen: neSII = {neSII}, TeOIII = {TeOIII}')

diags = pn.Diagnostics()
TeOIII, neSII = diags.getCrossTemDen(diag_tem='[OIII] 4363/5007',
                                     diag_den='[SII] 6731/6716',
                                     value_tem=df.loc['O3_4363A', 'int_epm']/df.loc['O3_5007A', 'int_epm'],
                                     value_den=df.loc['S2_6731A', 'int_epm']/df.loc['S2_6716A', 'int_epm'])
print(f'\nPyneb crossTemDen: neSII = {neSII}, TeOIII = {TeOIII}')


TeOII_Stasinska90 = (2 / ((1/(TeOIII/10000)) + 0.8)) * 10000.0
print(f'\nTeOII_Stasinska90 = {TeOII_Stasinska90}')


neSII = neSII
TeOII_Haegele06 = ((1.2 + 0.002*neSII + 4.2/neSII) / (10000.0/TeOIII + 0.08 + 0.003*neSII + 2.5/neSII)) * 10000.0
print(f'\nTeOII_Haegele06 = {TeOII_Haegele06}')


neSII = neSII
t_high = TeOIII/10000.0
TeOII_Haegele06 = ((1.2 + 0.002 * neSII + 4.2/neSII) / (1/t_high + 0.08 + 0.003*neSII + 2.5/neSII)) * 10000.0
print(f'\nTeOII_Haegele06 = {TeOII_Haegele06}')


ne = 412
T_low = TeOII_Haegele06
T_high = TeOIII
S2_abund = 12 + np.log10(S2.getIonAbundance(df.loc['S2_6716A', 'int_epm'], tem=T_low, den=ne, wave=6716, Hbeta=10000.))
O2_abund = 12 + np.log10(O2.getIonAbundance(df.loc['O2_3727A', 'int_epm'], tem=T_low, den=ne, to_eval='L(3726)+L(3729)', Hbeta=10000.))
N2_abund = 12 + np.log10(N2.getIonAbundance(df.loc['N2_6584A', 'int_epm'], tem=T_low, den=ne, wave=6584, Hbeta=10000.))
S3_abund = 12 + np.log10(S3.getIonAbundance(df.loc['S3_6312A', 'int_epm'], tem=T_low, den=ne, wave=6312, Hbeta=10000.))

O3_abund = 12 + np.log10(O3.getIonAbundance(df.loc['O3_5007A', 'int_epm'], tem=TeOIII, den=ne, wave=5007, Hbeta=10000.))

print('\n- Ionic abundances')
# print('S2_abund', S2_abund)
print('O2_abund', O2_abund)
# print('S3_abund', S3_abund)
print('O3_abund', O3_abund)
print('N2_abund', N2_abund)

T_low = 10900

O2_abund_epm = np.log10(df.loc['O2_3727A', 'int_epm']/10000) + 5.887 + 1.641/(T_low/10000) - 0.543*np.log10(T_low/10000) + 0.000114 * neSII
O3_abund_epm = np.log10((df.loc['O3_5007A', 'int_epm'] + df.loc['O3_4959A', 'int_epm'])/10000) + 6.1868 + 1.2491/(T_high/10000) - 0.5816*np.log10(T_high/10000)
N2_abund_epm = np.log10(((1 + 1/2.94) * df.loc['N2_6584A', 'int_epm'])/10000) + 6.291 + 0.90221/(T_low/10000) - 0.5511*np.log10(T_low/10000)

print('\nO2_abund_epm', O2_abund_epm)
print('O3_abund_epm', O3_abund_epm)
print('N2_abund_epm', N2_abund_epm)

# cm = Standard_DirectMetchod(n_steps=1000)
# flux_dict = cm.declare_line_fluxes(df.index.values, df['int_epm'].values, df['int_err'].values)
# int_dict = flux_dict
# cm.electron_diagnostics(int_dict, Thigh_diag=obsData[f'{obj}_chemical_model']['Te_high_diag'])
# cm.display_results()


# cHbeta, cHbeta_err = 0.017, 0.010
# red_law = obsData['sample_data']['red_law']
# RV = obsData['sample_data']['RV']
#
# rc = pn.RedCorr(R_V=RV, law=red_law, cHbeta=cHbeta)
# cors = rc.getCorr(df.wavelength.values, rel_wave=4861.0)
# lineInts = (df.flux.values / df.loc['H1_4861A'].flux) * cors
#
# df['int'] = lineInts * 10000
#
# # Plot diagnostics diagram for galaxy
# cm = Standard_DirectMetchod(n_steps=100)
# flux_dict = cm.declare_line_fluxes(df.index.values, df['flux'].values, df['eflux'].values)
#
# # Establish extinction correction
# f_Hbeta, f_lamdba = compute_arms_flambda(df, red_law, RV, ref_line='H1_4861A')
# int_dict = cm.red_corr(flux_dict, cHbeta=cHbeta, cHbeta_err=cHbeta_err, f_lambda=f_lamdba)
# cm.electron_diagnostics(int_dict, Thigh_diag=obsData[f'{obj}_chemical_model']['Te_high_diag'])
# cm.ionic_abundances(int_dict, obsData[f'{obj}_line_fitting'], obsData[f'{obj}_chemical_model'])
# cm.display_results()
#
# df['int'] = (df.flux.values / df.loc['H1_4861A'].flux) * np.power(10, cHbeta * f_lamdba)
# df['int_epm'] = df['int_epm'] * np.power(10, cHbeta * f_lamdba) / 10000
#
#
# R_O3 = (df.loc['O3_4959A'].int_epm + df.loc['O3_5007A'].int_epm) / df.loc['O3_4363A'].int_epm
# Te_O3_pm2017 = (0.7840 - 0.0001357 * R_O3 + 48.44 / R_O3) * 10000
# print('Te_O3_pm2017', Te_O3_pm2017)
#
# Te_S3_pm2017 = (1.19 * (Te_O3_pm2017 / 10000.0) - 0.32) * 10000
# print('Te_S3_pm2017', Te_S3_pm2017)
#
# R_S2 = df.loc['S2_6717A'].int_epm / df.loc['S2_6731A'].int_epm
# Te = 15618
# a0 = 16.054 - 7.79 / Te - 11.32 * Te
# a1 = -22.66 + 11.08 / Te + 16.02 * Te
# b0 = -21.61 + 11.89 / Te + 14.59 * Te
# b1 = 9.17 - 5.09 / Te - 6.18 * Te
# ne_S2_pm2017 = 1000 * (R_S2 * a0 + a1) / (R_S2 * b0 + b1)
# print('ne_S2_pm2017', ne_S2_pm2017)
#
# R_O2 = df.loc['S2_6717A'].int_epm / df.loc['S2_6731A'].int_epm
# ne = 412.0
# Te = 15618.0
# TO2_Hag2006 = ((1.2 + 0.002 * ne + 4.2/ne) / ((1/(Te/10000.0)) + 0.08 + 0.003*ne + 2.5/ne))*10000.0
# print('TO2_Hag2006', TO2_Hag2006)
#
# O3_coeff = (df.loc['O3_4959A', 'int_epm'] + df.loc['O3_5007A', 'int_epm']) / df.loc['H1_4861A', 'int_epm']
# O3_abund = np.log10(O3_coeff) + 6.1868 + 1.2491/(Te/10000.0) - 0.5816 * np.log10(Te/10000.0)
# print('O3_abund', O3_abund)
#
# T_low = TO2_Hag2006
# S2_coeff = (df.loc['S2_6717A', 'int_epm'] + df.loc['S2_6731A', 'int_epm']) / df.loc['H1_4861A', 'int_epm']
# S2_abund = np.log10(S2_coeff) + 5.463 + 0.941/(T_low/10000.0) - 0.97 * np.log10(T_low/10000.0)
# print('S2_abund', S2_abund)
