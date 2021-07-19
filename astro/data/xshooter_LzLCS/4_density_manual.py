from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import src.specsiser as sr
import pyneb as pn

c_KMpS = 299792.458

obsData = sr.loadConfData('./xshooter_LzLCS.ini')
data_folder = Path(obsData['data_location']['data_folder'])
results_folder = Path(obsData['data_location']['results_folder'])
objfile_list = obsData['data_location']['objfile_list']
sigmafile_list = obsData['data_location']['sigmafile_list']
objRef_list = obsData['data_location']['ref_list']
maskfile = obsData['data_location']['generalMask']

wmin_array = obsData['sample_data']['w_min_array']
wmax_array = obsData['sample_data']['w_max_array']
norm_flux = obsData['sample_data']['norm_flux']
z_obj = obsData['sample_data']['z_obj']
profile_conf = obsData['line_fitting']

DF_list = [None, None]
for i, objName in enumerate(objRef_list):

    # input data
    lineslog_file = results_folder/f'{objName}_linesLog.txt'

    # Load data
    linesDF = sr.lineslogFile_to_DF(lineslog_file)
    DF_list[i] = linesDF

# Join the dataframes
objDF = pd.concat(DF_list)

S2 = pn.Atom('S', 2)
O2 = pn.Atom('O', 2)

temp = 10000.0

ROII_narrow = objDF.loc['O2_3729A', 'gauss_flux'] / objDF.loc['O2_3726A', 'gauss_flux']
ROII_wide = objDF.loc['O2_3729A_w1', 'gauss_flux'] / objDF.loc['O2_3726A_w1', 'gauss_flux']

RSII_narrow = objDF.loc['S2_6716A', 'gauss_flux'] / objDF.loc['S2_6731A', 'gauss_flux']
RSII_wide = objDF.loc['S2_6716A_w1', 'gauss_flux'] / objDF.loc['S2_6731A_w1', 'gauss_flux']

neSII_narrow = S2.getTemDen(RSII_narrow, tem=temp, to_eval='L(6717)/L(6731)')
neSII_wide = S2.getTemDen(RSII_wide, tem=temp, to_eval='L(6717)/L(6731)')

neOII_narrow = O2.getTemDen(ROII_narrow, tem=temp, to_eval='L(3729)/L(3726)')
neOII_wide = O2.getTemDen(ROII_wide, tem=temp, to_eval='L(3729)/L(3726)')

print()
print('neSII_narrow', neSII_narrow)
print('neOII_narrow', neOII_narrow)
print()
print('neSII_wide', neSII_wide)
print('neOII_wide', neOII_wide)





# temp = 10000.0
# den = 200
#
# RSII_narrow = objDF.loc['S2_6716A', 'gauss_flux'] / objDF.loc['S2_6731A', 'gauss_flux']
# neSII_narrow = S2.getTemDen(RSII_narrow, tem=temp, to_eval='L(6717)/L(6731)')
# print(neSII_narrow)
#
# RSII_wide = objDF.loc['S2_6716A_w1', 'gauss_flux'] / objDF.loc['S2_6731A_w1', 'gauss_flux']
# neSII_wide = S2.getTemDen(RSII_wide, tem=temp, to_eval='L(6717)/L(6731)')
# print(neSII_wide)
# print(217 - 172)
# print('narrows difference', 172 - neSII_narrow)
# print('wides difference', 217 - neSII_wide)
#
# S2_emis_narrow = S2.getEmissivity(temp, den, wave=6716)/S2.getEmissivity(temp, den, wave=6731)
# S2_emis_wide = S2.getEmissivity(temp, den, wave=6716)/S2.getEmissivity(temp, den, wave=6731)
#
# O2_emis_narrow = O2.getEmissivity(temp, den, wave=3729)/O2.getEmissivity(temp, den, wave=3726)
# O2_emis_wide = O2.getEmissivity(temp, den, wave=3729)/O2.getEmissivity(temp, den, wave=3726)
#
# # S2_emis_narrow = S2.getEmissivity(temp, neSII_narrow, wave=6716)/S2.getEmissivity(temp, neSII_narrow, wave=6731)
# # S2_emis_wide = S2.getEmissivity(temp, neSII_wide, wave=6716)/S2.getEmissivity(temp, neSII_wide, wave=6731)
# #
# # O2_emis_narrow = O2.getEmissivity(temp, neSII_narrow, wave=3729)/O2.getEmissivity(temp, neSII_narrow, wave=3726)
# # O2_emis_wide = O2.getEmissivity(temp, neSII_wide, wave=3729)/O2.getEmissivity(temp, neSII_wide, wave=3726)
#
#
# print(f'Emissivity ratio [SII]6716/6731 narrow: {S2_emis_narrow:0.3f}')
# print(f'Emissivity ratio [SII]6716/6731 wide: {S2_emis_wide:0.3f}')
#
# print(f'Emissivity ratio [OII]3729/3726 narrow: {O2_emis_narrow:0.3f}')
# print(f'Emissivity ratio [OII]3729/3726 wide: {O2_emis_wide:0.3f}')