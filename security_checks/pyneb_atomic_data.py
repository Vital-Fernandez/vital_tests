import pyneb as pn

He1 = pn.RecAtom('He', 1)
O3 = pn.Atom('O', 3)

print(pn.__version__)
He1.printSources()

print('O3 line ratio', O3.getEmissivity(8000, 1000, wave=5007) / O3.getEmissivity(8000, 1000, wave=4959))


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

