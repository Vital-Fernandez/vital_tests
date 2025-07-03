import pyneb as pn
from matplotlib import pyplot as plt
He1 = pn.RecAtom('He', 1)
H1 = pn.RecAtom('H', 1)
O3 = pn.Atom('O', 3)
Ar4 = pn.Atom('Ar', 4)
print(He1.printTransition(5876))


print('Halpha/Hbeta', H1.getEmissivity(10000, 100, wave=6563)/H1.getEmissivity(10000, 100, wave=4861))


S3 = pn.Atom('S', 3)
H1.getWave(14, 2)
H1.getWave(13, 2)
N2 = pn.Atom('N', 2)
# N2.plotGrotrian()
plt.show()
print('S3 line ratio', S3.getEmissivity(8000, 1000, wave=9531) / S3.getEmissivity(8000, 1000, wave=9069))
print('N2 line ratio', N2.getEmissivity(8000, 1000, wave=6583) / N2.getEmissivity(8000, 1000, wave=6548))
print(N2.printSources())
print('----------')
#
# print(pn.__version__)
# He1.printSources()
#
print('O3 line ratio', O3.getEmissivity(8000, 1000, wave=5007) / O3.getEmissivity(8000, 1000, wave=4959))
print(O3.printSources())
print('----------')
O3.plotGrotrian()
plt.show()
# print('Ar4 line ratio', Ar4.getEmissivity(15000, 100, wave=4740) / Ar4.getEmissivity(15000, 100, wave=4710))
#

# H1 = pn.RecAtom('H', 1)
# print('HBeta/Halpha Paschen', H1.getEmissivity(8000, 1000, wave=4861) / H1.getEmissivity(8000, 1000, wave=18751))
# print('HBeta/Halpha Paschen', H1.getEmissivity(20000, 10, wave=4861) / H1.getEmissivity(20000, 10, wave=18751))
# print('HBeta/Halpha Paschen', H1.getEmissivity(10000, 100, wave=4861) / H1.getEmissivity(10000, 100, wave=18751))

# coso = pn.getAtomDict(['Ar3', 'Ar4', 'He1'])
#
# print(coso)
#
# coso2 = pn.getAtomDict(['Ar3', 'Ar4', 'He1'], only_coll=True)
#
# print(coso2)

