import numpy as np
import pyneb as pn

H1 = pn.RecAtom('H', 1)

ne, Te = 100, 10000.0

print('Halpha/Hbeta', H1.getEmissivity(Te, ne, wave=6563) / H1.getEmissivity(Te, ne, wave=4861))
print('Halpha/Hgamma', H1.getEmissivity(Te, ne, wave=6563) / H1.getEmissivity(Te, ne, wave=4341))
print('Halpha/Hdelta', H1.getEmissivity(Te, ne, wave=6563) / H1.getEmissivity(Te, ne, wave=4101))

abund=0.005
abund2=0.005

O2 = pn.Atom('O', 2)
S2 = pn.Atom('S', 2)
O3 = pn.Atom('O', 3)
S3 = pn.Atom('S', 3)
N2 = pn.Atom('N', 2)


Hbeta_em = H1.getEmissivity(Te, ne, wave=4861)
S2_6717A_em = S2.getEmissivity(Te, ne, wave=6717)
S2_6731_em = S2.getEmissivity(Te, ne, wave=6731)

# print(S2_6717A_em/S2_6731_em)

S2_6717A_flux = abund * S2_6717A_em/Hbeta_em
S2_6731_flux = abund * S2_6731_em/Hbeta_em

S3_6312A_em = S3.getEmissivity(Te, ne, wave=6312)
S3_9069A_em = S3.getEmissivity(Te, ne, wave=9069)
S3_9531A_em = S3.getEmissivity(Te, ne, wave=9531)

S3_6312A_flux = abund2 * S3_6312A_em/Hbeta_em
S3_9069A_flux = abund2 * S3_9069A_em/Hbeta_em
S3_9531A_flux = abund2 * S3_9531A_em/Hbeta_em

S2_4069A_em = S2.getEmissivity(Te, ne, wave=4069)
S2_4076A_em = S2.getEmissivity(Te, ne, wave=4076)

# print(S2_4069A_em/S2_4076A_em)
# print(S2.getEmissivity(18500, ne, wave=4069) / S2.getEmissivity(18500, ne, wave=4076))
# print(S2.getEmissivity(18500, 500, wave=4069) / S2.getEmissivity(18500, 500, wave=4076))
# print(S2.getEmissivity(Te, 500, wave=4069) / S2.getEmissivity(Te, 500, wave=4076))
# print(S2.getEmissivity(8000, ne, wave=4069) / S2.getEmissivity(8000, ne, wave=4076))

# S2_6717A_array = np.random.normal(S2_6717A_mean, S2_6717A_err, size=1000)
# S2_6731A_array = np.random.normal(S2_6731A_mean, S2_6731A_err, size=1000)
# RSII_array = S2_6717A_array/S2_6731A_array

# #NARROW, BROAD, GLOBAL
# Sii17_j1152 = np.array([15.61,12.11,27.73])
# Sii17_j1152_err = np.array([1.84,2.869,3.409])
# Sii31_j1152 = np.array([13.59,14.13,27.71])
# Sii31_j1152_err = np.array([1.454,2.235,2.667])

n_iteration = 1000

S2_6717_n1 = np.random.normal(15.61, 1.84, size=n_iteration)
S2_6717_B = np.random.normal(12.11, 2.869, size=n_iteration)
S2_6717_G = np.random.normal(27.73, 3.409, size=n_iteration)

S2_6731_n1 = np.random.normal(13.59, 1.454, size=n_iteration)
S2_6731_B = np.random.normal(14.13, 2.235, size=n_iteration)
S2_6731_G = np.random.normal(27.71, 2.667, size=n_iteration)

temp_J1152 = np.random.normal(13430.0, 900, size=n_iteration)

nSII_n1 = S2.getTemDen(S2_6717_n1/S2_6731_n1, tem=temp_J1152, to_eval='L(6717)/L(6731)')
nSII_B = S2.getTemDen(S2_6717_B/S2_6731_B, tem=temp_J1152, to_eval='L(6717)/L(6731)')
nSII_G = S2.getTemDen(S2_6717_G/S2_6731_G, tem=temp_J1152, to_eval='L(6717)/L(6731)')

print(np.nanmean(nSII_n1), np.nanstd(nSII_n1))
print(np.nanmean(nSII_B), np.nanstd(nSII_B), np.isnan(nSII_B).sum())
print(np.nanmean(nSII_G), np.nanstd(nSII_G))

ne_SII_ = S2.getTemDen(np.array([15.61,12.11,27.73]) / np.array([13.59,14.13,27.71]), tem=13430.0, to_eval='L(6717)/L(6731)')

print(ne_SII_)

print()
print()


S2_6717_n1 = np.random.normal(16.73, 0.827, size=n_iteration)
S2_6717_n2 = np.random.normal(16.07, 0.639, size=n_iteration)
S2_6717_B = np.random.normal(3.74, 1.591, size=n_iteration)
S2_6717_G = np.random.normal(36.54, 1.903, size=n_iteration)

# Sii31_j0925 = np.array([13.81,13.42,3.93,31.16])
# Sii31_j0925_err = np.array([1.005,0.860,2.029,2.422])

S2_6731_n1 = np.random.normal(13.81, 1.005, size=n_iteration)
S2_6731_n2 = np.random.normal(13.42, 0.860, size=n_iteration)
S2_6731_B = np.random.normal(3.93, 2.029, size=n_iteration)
S2_6731_G = np.random.normal(31.16, 2.422, size=n_iteration)

temp_J0925 = np.random.normal(15010.0, 410, size=n_iteration)
temp_J0925_B = np.random.normal(10000.0, 500, size=n_iteration)

nSII_n1 = S2.getTemDen(S2_6717_n1/S2_6731_n1, tem=temp_J0925, to_eval='L(6717)/L(6731)')
nSII_n2 = S2.getTemDen(S2_6717_n2/S2_6731_n2, tem=temp_J0925, to_eval='L(6717)/L(6731)')
nSII_B = S2.getTemDen(S2_6717_B/S2_6731_B, tem=temp_J0925_B, to_eval='L(6717)/L(6731)')
nSII_G = S2.getTemDen(S2_6717_G/S2_6731_G, tem=temp_J0925, to_eval='L(6717)/L(6731)')


print(np.nanmean(nSII_n1), np.nanstd(nSII_n1))
print(np.nanmean(nSII_n2), np.nanstd(nSII_n2))
print(np.nanmean(nSII_B), np.nanstd(nSII_B), np.isnan(nSII_B).sum())
print(np.nanmean(nSII_G), np.nanstd(nSII_G))

ne_SII_ = S2.getTemDen(np.array([16.73, 16.07, 3.74, 36.54]) / np.array([13.81, 13.42, 3.93, 31.16]), tem=15010.0, to_eval='L(6717)/L(6731)')
print(ne_SII_)

# ne_SII = S2.getTemDen(RSII_array, tem=10000.0, to_eval='L(6717)/L(6731)')
# ne_SII_inv = S2.getTemDen(S2_6731_flux/S2_6717A_flux, tem=11000.0, to_eval='L(6731)/L(6717)')
# print(ne_SII)
# print(ne_SII_inv, '\n')

# print('This flux example ', S2_6717A_flux/S2_6731_flux)
# print('This temperature', S2.getTemDen(int_ratio=1.50, tem=10000.0, to_eval='L(6716)/L(6731)'))
#
#
# S2_EG = pn.EmisGrid('S', 2,  n_tem=100, n_den=100,
#                     tem_min=5000., tem_max=50000., den_min=0.01,
#                     den_max=1000, restore_file=None, atomObj=None)
# S2_EG.plotContours(to_eval = 'L(6716)/L(6731)', log_levels=False)
# plt.show()

# plt.show()
diags = pn.Diagnostics()
diags.addDiag("[SIII] 6312/9531", diag_tuple=('S3', 'L(6312)/L(9531)', 'RMS([E(9069), E(9531)])'))
diags.addDiag("[SII] 6716/6731", diag_tuple=('S2', 'L(6716)/L(6731)', 'RMS([E(6716), E(6731)])'))
#
#
# attr1 = dict(diag_den = "[SII] 6716/6731",
#              value_den = S2_6717A_flux/S2_6731_flux,
#              diag_tem = '[SIII] 6312/9200+',
#              value_tem = S3_6312A_flux/(S3_9069A_flux+S3_9531A_flux))
#
# attr2 = dict(diag_den = '[SII] 6731/6716',
#              value_den = S2_6731_flux/S2_6717A_flux,
#              diag_tem = '[SIII] 6312/9531',
#              value_tem = S3_6312A_flux/S3_9531A_flux)
#
# attr3 = dict(diag_den = '[SII] 6731/6716',
#              value_den = S2_6731_flux/S2_6717A_flux,
#              diag_tem = '[SIII] 6312/9069',
#              value_tem = S3_6312A_flux/S3_9069A_flux)
#
# # diags2 = pn.Diagnostics()
# # diags2.addDiag("[SII] 6731/6716", diag_tuple=('S2', 'L(6731)/L(6716)', 'RMS([E(6731), E(6716)])'))
# # diags2.addDiag("[SIII] 6312/9531", diag_tuple=('S3', 'L(6312)/L(9531)', 'RMS([E(9069), E(9531)])'))
#
# tem_array, den_array = diags.getCrossTemDen(**attr1)
# tem_array2, den_array2 = diags.getCrossTemDen(**attr2)
# tem_array3, den_array3 = diags.getCrossTemDen(**attr3)
#
# print(tem_array, den_array)
# print(tem_array2, den_array2)
# print(tem_array3, den_array3)




# import re
#
# print(re.findall(pattern=r'\d+', string="L(6312)/L(9531)"))
# print(re.findall(pattern=r'\d+', string='L(4626)/(L(6435)+L(7005))'))
#

# Te_value = 5000
# ne_SII_ratio = S2.getEmissivity(tem=Te_value, den=ne_lowReg, wave=6716)/S2.getEmissivity(tem=Te_value, den=ne_lowReg, wave=6731)
# den_SII_dist = S2.getTemDen(ne_SII_ratio, tem=10000.0, to_eval='L(6716)/L(6731)')
# print(den_SII_dist.mean(), den_SII_dist.std())
# diags = pn.Diagnostics()



# print(diags.getCrossTemDen(**attr2))






# line_labels = ['S2_6716A', 'S2_6731A', 'O3_4363A', 'O3_4959A', 'O3_5007A']
# obs = pn.Observation()
# for lineLabel in line_labels:
#     newLine = pn.EmissionLine(label=lineLabel, obsIntens=np.nan)
#     obs.addLine(newLine)
# diags = pn.Diagnostics()
# diags.addDiagsFromObs(obs)
# print(diags.getDiags())
# diags.getCrossTemDen()

# Hbeta_em = H1.getEmissivity(Te, ne, wave=4861)
# O2_3726A_em = O2.getEmissivity(Te, ne, wave=3726)
# O2_3729A_em = O2.getEmissivity(Te, ne, wave=3729)
# O2_3720b_em = O2_3726A_em + O2_3729A_em
#
# O2_3726A_flux = abund * O2_3726A_em/Hbeta_em
# O2_3729A_flux = abund * O2_3729A_em/Hbeta_em
# O2_3720b_flux = abund * O2_3720b_em/Hbeta_em


# S2.printIonic()
#
# Hbeta_em = H1.getEmissivity(Te, ne, wave=4861)
# O2_3726A_em = O2.getEmissivity(Te, ne, wave=3726)
# O2_3729A_em = O2.getEmissivity(Te, ne, wave=3729)
# O2_3720b_em = O2_3726A_em + O2_3729A_em
#
# O2_3726A_flux = abund * O2_3726A_em/Hbeta_em
# O2_3729A_flux = abund * O2_3729A_em/Hbeta_em
# O2_3720b_flux = abund * O2_3720b_em/Hbeta_em
#
# O2_3726A_ref = 'I(3,1)'
# O2_3729A_ref = 'I(2,1)'
# O2_3720b_ref = 'I(2,1)+I(3,1)'
#
# O2_3726A_abund = O2.getIonAbundance(O2_3726A_flux, Te, ne, to_eval=O2_3726A_ref, Hbeta=1)
# O2_3729A_abund = O2.getIonAbundance(O2_3729A_flux, Te, ne, to_eval=O2_3729A_ref, Hbeta=1)
# O2_3720b_abund = O2.getIonAbundance(O2_3720b_flux, Te, ne, to_eval=O2_3720b_ref, Hbeta=1)
#
# print('O2_3726A_abund', O2_3726A_abund)
# print('O2_3729A_abund', O2_3729A_abund)
# print('O2_3720b_abund', O2_3729A_abund)



# "OII 4649.13/4089.29" : ('O2r', "S('4649.13')/S('4089.29')", "RMS([SE('4649.13'), SE('4089.29')])")
# "OII 4649.13/4661.63" : ('O2r', "S('4649.13')/S('4661.63')", "RMS([SE('4649.13'), SE('4661.63')])")
# "[ArIII] (7751+7136)/9m" : ('Ar3', '(L(7751)+L(7136))/L(90000)', 'RMS([E(90000), E(7751)*L(7751)/(L(7751)+L(7136)), E(7136)*L(7136)/(L(7751)+L(7136))])')
# "[ArIII] 5192/7136" : ('Ar3', 'L(5192)/L(7136)', 'RMS([E(7136), E(5192)])')
# "[ArIII] 5192/7300+" : ('Ar3', 'L(5192)/(L(7751)+L(7136))', 'RMS([E(7751)*L(7751)/(L(7751)+L(7136)), E(7136)*L(7136)/(L(7751)+L(7136)), E(5192)])')
# "[ArIII] 7136/9m" : ('Ar3', 'L(7136)/L(90000)', 'RMS([E(90000), E(7136)])')
# "[ArIII] 9.0m/21.8m" : ('Ar3', 'L(89897)/L(218000)', 'RMS([E(89897), E(218000)])')
# "[ArIV] 2860+/4720+" : ('Ar4', '(L(2854)+L(2868))/(L(4711)+L(4740))', 'RMS([E(4711)*L(4711)/(L(4711)+L(4740)), E(4740)*L(4740)/(L(4711)+L(4740)), E(2854)*L(2854)/(L(2854)+L(2868)), E(2868)*L(2854)/(L(2854)+L(2868))])')
# "[ArIV] 4740/4711" : ('Ar4', 'L(4740)/L(4711)', 'RMS([E(4711), E(4740)])')
# "[ArIV] 7230+/4720+" : ('Ar4', '(L(7170)+L(7263))/(L(4711)+L(4740))', 'RMS([E(4711)*L(4711)/(L(4711)+L(4740)), E(4740)*L(4740)/(L(4711)+L(4740)), E(7170)*L(7170)/(L(7170)+L(7263)), E(7263)*L(7263)/(L(7170)+L(7263))])')
# "[ArV] 4626/6600+" : ('Ar5', 'L(4626)/(L(6435)+L(7005))', 'RMS([E(6435)*L(6435)/(L(6435)+L(7005)), E(7005)*L(7005)/(L(6435)+L(7005)), E(4626)])')
# "[CIII] 1909/1907" : ('C3', 'L(1909)/L(1907)', 'RMS([E(1909), E(1907)])')
# "[ClIII] 5538/5518" : ('Cl3', 'L(5538)/L(5518)', 'RMS([E(5518), E(5538)])')
# "[ClIV] 5323/7531" : ('Cl4', 'L(5323)/L(7531)', 'RMS([E(7531), E(5323)])')
# "[ClIV] 5323/7700+" : ('Cl4', 'L(5323)/(L(7531)+L(8046))', 'RMS([E(7531)*L(7531)/(L(7531)+L(8046)), E(8046)*L(8046)/(L(7531)+L(8046)), E(5323)])')
# "[FeIII] 4659/4009" : ('Fe3', 'L(4659)/L(4009)', 'RMS([E(4659), E(4009)])')
# "[FeIII] 4659/4701" : ('Fe3', 'L(4659)/L(4701)', 'RMS([E(4659), E(4701)])')
# "[FeIII] 4659/4734" : ('Fe3', 'L(4659)/L(4734)', 'RMS([E(4659), E(4734)])')
# "[FeIII] 4701/4009" : ('Fe3', 'L(4701)/L(4009)', 'RMS([E(4701), E(4009)])')
# "[FeIII] 4701/4734" : ('Fe3', 'L(4701)/L(4734)', 'RMS([E(4701), E(4734)])')
# "[FeIII] 4734/4009" : ('Fe3', 'L(4734)/L(4009)', 'RMS([E(4734), E(4009)])')
# "[FeIII] 4881/4009" : ('Fe3', 'L(4881)/L(4009)', 'RMS([E(4881), E(4009)])')
# "[FeIII] 4881/4659" : ('Fe3', 'L(4881)/L(4659)', 'RMS([E(4881), E(4659)])')
# "[FeIII] 4881/4701" : ('Fe3', 'L(4881)/L(4701)', 'RMS([E(4881), E(4701)])')
# "[FeIII] 4881/4734" : ('Fe3', 'L(4881)/L(4734)', 'RMS([E(4881), E(4734)])')
# "[FeIII] 4881/4931" : ('Fe3', 'L(4881)/L(4931)', 'RMS([E(4881), E(4931)])')
# "[FeIII] 4881/5011" : ('Fe3', 'L(4881)/L(5011)', 'RMS([E(4881), E(5011)])')
# "[FeIII] 4925/4009" : ('Fe3', 'L(4925)/L(4009)', 'RMS([E(4925), E(4009)])')
# "[FeIII] 4925/4659" : ('Fe3', 'L(4925)/L(4659)', 'RMS([E(4925), E(4659)])')
# "[FeIII] 4925/4701" : ('Fe3', 'L(4925)/L(4701)', 'RMS([E(4925), E(4701)])')
# "[FeIII] 4925/4734" : ('Fe3', 'L(4925)/L(4734)', 'RMS([E(4925), E(4734)])')
# "[FeIII] 4925/4881" : ('Fe3', 'L(4925)/L(4881)', 'RMS([E(4925), E(4881)])')
# "[FeIII] 4925/4931" : ('Fe3', 'L(4925)/L(4931)', 'RMS([E(4925), E(4931)])')
# "[FeIII] 4925/5011" : ('Fe3', 'L(4925)/L(5011)', 'RMS([E(4925), E(5011)])')
# "[FeIII] 4931/4009" : ('Fe3', 'L(4931)/L(4009)', 'RMS([E(4931), E(4009)])')
# "[FeIII] 4931/4659" : ('Fe3', 'L(4931)/L(4659)', 'RMS([E(4931), E(4659)])')
# "[FeIII] 4931/4701" : ('Fe3', 'L(4931)/L(4701)', 'RMS([E(4931), E(4701)])')
# "[FeIII] 4931/4734" : ('Fe3', 'L(4931)/L(4734)', 'RMS([E(4931), E(4734)])')
# "[FeIII] 4987/4009" : ('Fe3', 'L(4987)/L(4009)', 'RMS([E(4987), E(4009)])')
# "[FeIII] 4987/4659" : ('Fe3', 'L(4987)/L(4659)', 'RMS([E(4987), E(4659)])')
# "[FeIII] 4987/4701" : ('Fe3', 'L(4987)/L(4701)', 'RMS([E(4987), E(4701)])')
# "[FeIII] 4987/4734" : ('Fe3', 'L(4987)/L(4734)', 'RMS([E(4987), E(4734)])')
# "[FeIII] 4987/4881" : ('Fe3', 'L(4987)/L(4881)', 'RMS([E(4987), E(4881)])')
# "[FeIII] 4987/4925" : ('Fe3', 'L(4987)/L(4925)', 'RMS([E(4987), E(4925)])')
# "[FeIII] 4987/4931" : ('Fe3', 'L(4987)/L(4931)', 'RMS([E(4987), E(4931)])')
# "[FeIII] 4987/5011" : ('Fe3', 'L(4987)/L(5011)', 'RMS([E(4987), E(5011)])')
# "[FeIII] 5011/4009" : ('Fe3', 'L(5011)/L(4009)', 'RMS([E(5011), E(4009)])')
# "[FeIII] 5011/4659" : ('Fe3', 'L(5011)/L(4659)', 'RMS([E(5011), E(4659)])')
# "[FeIII] 5011/4701" : ('Fe3', 'L(5011)/L(4701)', 'RMS([E(5011), E(4701)])')
# "[FeIII] 5011/4734" : ('Fe3', 'L(5011)/L(4734)', 'RMS([E(5011), E(4734)])')
# "[FeIII] 5011/4931" : ('Fe3', 'L(5011)/L(4931)', 'RMS([E(5011), E(4931)])')
# "[FeIII] 5270/4009" : ('Fe3', 'L(5270)/L(4009)', 'RMS([E(5270), E(4009)])')
# "[FeIII] 5270/4659" : ('Fe3', 'L(5270)/L(4659)', 'RMS([E(5270), E(4659)])')
# "[FeIII] 5270/4701" : ('Fe3', 'L(5270)/L(4701)', 'RMS([E(5270), E(4701)])')
# "[FeIII] 5270/4734" : ('Fe3', 'L(5270)/L(4734)', 'RMS([E(5270), E(4734)])')
# "[FeIII] 5270/4881" : ('Fe3', 'L(5270)/L(4881)', 'RMS([E(5270), E(4881)])')
# "[FeIII] 5270/4925" : ('Fe3', 'L(5270)/L(4925)', 'RMS([E(5270), E(4925)])')
# "[FeIII] 5270/4931" : ('Fe3', 'L(5270)/L(4931)', 'RMS([E(5270), E(4931)])')
# "[FeIII] 5270/4987" : ('Fe3', 'L(5270)/L(4987)', 'RMS([E(5270), E(4987)])')
# "[FeIII] 5270/5011" : ('Fe3', 'L(5270)/L(5011)', 'RMS([E(5270), E(5011)])')
# "[NIII] 1751+/57.4m" : ('N3', 'B("1751A+")/L(574000)', 'RMS([E(574000), BE("1751A+")])')
# "[NII] 121m/20.5m" : ('N2', 'L(1214747)/L(2054427)', 'RMS([E(2054427)/E(1214747)])')
# "[NII] 5755/6548" : ('N2', 'L(5755)/L(6548)', 'RMS([E(6548), E(5755)])')
# "[NII] 5755/6584" : ('N2', 'L(5755)/L(6584)', 'RMS([E(6584), E(5755)])')
# "[NII] 5755/6584+" : ('N2', 'L(5755)/(L(6548)+L(6584))', 'RMS([E(6548)*L(6548)/(L(6548)+L(6584)), E(6584)*L(6584)/(L(6584)+L(6548)), E(5755)])')
# "[NI] 5198/5200" : ('N1', 'I(3, 1)/I(2, 1)', 'RMS([E(5200), E(5198)])')
# "[NeIII] 15.6m/36.0m" : ('Ne3', 'L(156000)/L(360000)', 'RMS([E(156000), E(360000)])')
# "[NeIII] 3343/3930+" : ('Ne3', 'L(3343)/(L(3869)+L(3968))', 'RMS([E(3869)*L(3869)/(L(3869)+L(3968)), E(3968)*L(3968)/(L(3869)+L(3968)), E(3343)])')
# "[NeIII] 3344/3930+" : ('Ne3', 'L(3343)/(L(3869)+L(3968))', 'RMS([E(3869)*L(3869)/(L(3869)+L(3968)), E(3968)*L(3968)/(L(3869)+L(3968)), E(3343)])')
# "[NeIII] 3869/15.6m" : ('Ne3', 'L(3869)/L(156000)', 'RMS([E(156000), E(3869)])')
# "[NeIII] 3930+/15.6m" : ('Ne3', '(L(3869)+L(3968))/L(156000)', 'RMS([E(156000), E(3869)*L(3869)/(L(3869)+L(3968)), E(3968)*L(3968)/(L(3869)+L(3968))])')
# "[NeV] 14.3m/24.2m" : ('Ne5', 'L(143000)/L(242000)', 'RMS([E(143000), E(242000)])')
# "[NeV] 1575/3426" : ('Ne5', 'L(1575)/L(3426)', 'RMS([E(1575), E(3426)])')
# "[NeV] 2973/3370+" : ('Ne5', 'L(2973)/(L(3426)+L(3346))', 'RMS([E(3426)*L(3426)/(L(3426)+L(3346)), E(3346)*L(3346)/(L(3426)+L(3346)), E(2973)])')
# "[NeV] 3426/24.2m" : ('Ne5', 'L(3426)/L(242000)', 'RMS([E(3426), E(242000)])')
# "[NiIII] 6000/6401" : ('Ni3', 'L(6000)/L(6401)', 'RMS([E(6000), E(6401)])')
# "[NiIII] 6000/6682" : ('Ni3', 'L(6000)/L(6682)', 'RMS([E(6000), E(6682)])')
# "[NiIII] 6000/6797" : ('Ni3', 'L(6000)/L(6797)', 'RMS([E(6000), E(6797)])')
# "[NiIII] 6000/7125" : ('Ni3', 'L(6000)/L(7125)', 'RMS([E(6000), E(7125)])')
# "[NiIII] 6000/7890" : ('Ni3', 'L(6000)/L(7890)', 'RMS([E(6000), E(7890)])')
# "[NiIII] 6000/8500" : ('Ni3', 'L(6000)/L(8500)', 'RMS([E(6000), E(8500)])')
# "[NiIII] 6401/7125" : ('Ni3', 'L(6401)/L(7125)', 'RMS([E(6401), E(7125)])')
# "[NiIII] 6401/7890" : ('Ni3', 'L(6401)/L(7890)', 'RMS([E(6401), E(7890)])')
# "[NiIII] 6401/8500" : ('Ni3', 'L(6401)/L(8500)', 'RMS([E(6401), E(8500)])')
# "[NiIII] 6534/6401" : ('Ni3', 'L(6534)/L(6401)', 'RMS([E(6534), E(6401)])')
# "[NiIII] 6534/6682" : ('Ni3', 'L(6534)/L(6682)', 'RMS([E(6534), E(6682)])')
# "[NiIII] 6534/6797" : ('Ni3', 'L(6534)/L(6797)', 'RMS([E(6534), E(6797)])')
# "[NiIII] 6534/7125" : ('Ni3', 'L(6534)/L(7125)', 'RMS([E(6534), E(7125)])')
# "[NiIII] 6534/7890" : ('Ni3', 'L(6534)/L(7890)', 'RMS([E(6534), E(7890)])')
# "[NiIII] 6534/8500" : ('Ni3', 'L(6534)/L(8500)', 'RMS([E(6534), E(8500)])')
# "[NiIII] 6682/7125" : ('Ni3', 'L(6682)/L(7125)', 'RMS([E(6682), E(7125)])')
# "[NiIII] 6682/7890" : ('Ni3', 'L(6682)/L(7890)', 'RMS([E(6682), E(7890)])')
# "[NiIII] 6682/8500" : ('Ni3', 'L(6682)/L(8500)', 'RMS([E(6682), E(8500)])')
# "[NiIII] 6797/7125" : ('Ni3', 'L(6797)/L(7125)', 'RMS([E(6797), E(7125)])')
# "[NiIII] 6797/7890" : ('Ni3', 'L(6797)/L(7890)', 'RMS([E(6797), E(7890)])')
# "[NiIII] 6797/8500" : ('Ni3', 'L(6797)/L(8500)', 'RMS([E(6797), E(8500)])')
# "[NiIII] 6946/6401" : ('Ni3', 'L(6946)/L(6401)', 'RMS([E(6946), E(6401)])')
# "[NiIII] 6946/6682" : ('Ni3', 'L(6946)/L(6682)', 'RMS([E(6946), E(6682)])')
# "[NiIII] 6946/6797" : ('Ni3', 'L(6946)/L(6797)', 'RMS([E(6946), E(6797)])')
# "[NiIII] 6946/7125" : ('Ni3', 'L(6946)/L(7125)', 'RMS([E(6946), E(7125)])')
# "[NiIII] 6946/7890" : ('Ni3', 'L(6946)/L(7890)', 'RMS([E(6946), E(7890)])')
# "[NiIII] 6946/8500" : ('Ni3', 'L(6946)/L(8500)', 'RMS([E(6946), E(8500)])')
# "[OIII] 1664+/5007" : ('O3', '(B("1664A+"))/L(5007)', 'RMS([BE("1664A+"), E(5007)])')
# "[OIII] 1666/4363" : ('O3', 'L(1666)/L(4363)', 'RMS([E(4363), E(1666)])')
# "[OIII] 1666/5007" : ('O3', 'L(1666)/L(5007)', 'RMS([E(5007), E(1666)])')
# "[OIII] 1666/5007+" : ('O3', 'L(1666)/(L(5007)+L(4959))', 'RMS([E(5007)*L(5007)/(L(5007)+L(4959)), E(4959)*L(4959)/(L(5007)+L(4959)), E(1666)])')
# "[OIII] 4363/5007" : ('O3', 'L(4363)/L(5007)', 'RMS([E(5007), E(4363)])')
# "[OIII] 4363/5007+" : ('O3', 'L(4363)/(L(5007)+L(4959))', 'RMS([E(5007)*L(5007)/(L(5007)+L(4959)), E(4959)*L(4959)/(L(5007)+L(4959)), E(4363)])')
# "[OIII] 5007/88m" : ('O3', 'L(5007)/L(883000)', 'RMS([E(883000), E(5007)])')
# "[OIII] 51m/88m" : ('O3', 'L(518000)/L(883000)', 'RMS([E(883000), E(518000)])')
# "[OII] 3726/3729" : ('O2', 'L(3726)/L(3729)', 'RMS([E(3729), E(3726)])')
# "[OII] 3727+/7325+" : ('O2', '(L(3726)+L(3729))/(B("7319A+")+B("7330A+"))', 'RMS([E(3726)*L(3726)/(L(3726)+L(3729)), E(3729)*L(3729)/(L(3726)+L(3729)),BE("7319A+")*B("7319A+")/(B("7319A+")+B("7330A+")),BE("7330A+")*B("7330A+")/(B("7319A+")+B("7330A+"))])')
# "[OII] 3727+/7325+b" : ('O2', '(L(3726)+L(3729))/(I(4,2)+I(4,3)+I(5,2)+I(5,3))', 'RMS([E(3726)*L(3726)/(L(3726)+L(3729)), E(3729)*L(3729)/(L(3726)+L(3729)),BE("7319A+")*B("7319A+")/(B("7319A+")+B("7330A+")),BE("7330A+")*B("7330A+")/(B("7319A+")+B("7330A+"))])')
# "[OII] 3727+/7325+c" : ('O2', '(B("3727A+"))/(B("7319A+")+B("7330A+"))', 'RMS([BE("7319A+")*B("7319A+")/(B("7319A+")+B("7330A+")), BE("7330A+")*B("7330A+")/(B("7319A+")+B("7330A+")), BE("3727A+")])')
# "[OIV] 1400+/25.9m" : ('O4', 'B("1400A+")/L(259000)', 'RMS([BE("1400A+"), E(259000)])')
# "[OIV] 1401/1405" : ('O4', 'L(1401)/L(1405)', 'RMS([E(1401), E(1405)])')
# "[OI] 5577/6300" : ('O1', 'L(5577)/L(6300)', 'RMS([E(6300), E(5577)])')
# "[OI] 5577/6300+" : ('O1', 'L(5577)/(L(6300)+L(6364))', 'RMS([E(6300)*L(6300)/(L(6300)+L(6364)), E(6364)*L(6364)/(L(6300)+L(6364)), E(5577)])')
# "[OI] 5577/6302" : ('O1', 'L(5577)/L(6300)', 'RMS([E(6300), E(5577)])')
# "[OI] 63m/147m" : ('O1', 'L(632000)/L(1455000)', 'RMS([E(632000), E(1455000)])')
# "[SIII] 18.7m/33.5m" : ('S3', 'L(187000)/L(335000)', 'RMS([E(335000), E(187000)])')
# "[SIII] 6312/18.7m" : ('S3', 'L(6312)/L(187000)', 'RMS([E(187000), E(6312)])')
# "[SIII] 6312/9069" : ('S3', 'L(6312)/L(9069)', 'RMS([E(9069), E(6312)])')
# "[SIII] 6312/9200+" : ('S3', 'L(6312)/(L(9069)+L(9531))', 'RMS([E(9069)*L(9069)/(L(9069)+L(9531)), E(9531)*L(9531)/(L(9069)+L(9531)), E(6312)])')
# "[SIII] 9069/18.7m" : ('S3', 'L(9069)/L(187000)', 'RMS([E(187000), E(9069)])')
# "[SII] 4069/4076" : ('S2', 'L(4069)/L(4076)', 'RMS([E(4069), E(4076)])')
# "[SII] 4072+/6720+" : ('S2', '(L(4069)+L(4076))/(L(6716)+L(6731))', 'RMS([E(6716)*L(6716)/(L(6716)+L(6731)), E(6731)*L(6731)/(L(6716)+L(6731)), E(4069)*L(4069)/(L(4069)+L(4076)), E(4076)*L(4076)/(L(4069)+L(4076))])')
# "[SII] 6731/6716" : ('S2', 'L(6731)/L(6716)', 'RMS([E(6716), E(6731)])')
