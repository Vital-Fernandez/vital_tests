import matplotlib.pyplot as plt
import pyneb as pn

O3_EG = pn.EmisGrid('O', 3, n_tem=30, n_den=30)
O3_EG = pn.EmisGrid(elem='O', spec=3, n_tem=100, n_den=100,
                    tem_min=5000., tem_max=30000., den_min=10.,
                    den_max=1.e10, restore_file=None, atomObj=None)
O3_5007 = O3_EG.getGrid(wave=5007)
O3_Te = O3_EG.getGrid(to_eval = 'L(4363)/L(5007)')
f, ax = plt.subplots(figsize=(7,5))
O3_EG.plotContours(to_eval = 'L(4363)/L(5007)', ax=ax)
f.savefig('OIII_diag.png')