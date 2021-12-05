import pyneb as pn

Ar3 = pn.Atom('Ar', 3)
Ar4 = pn.Atom('Ar', 4)
He1 = pn.RecAtom('He', 1)
O2 = pn.Atom('O', 3)

O2.printSources()
# Ar3.printSources()
# Ar4.printSources()

# coso = pn.getAtomDict(['Ar3', 'Ar4', 'He1'])
#
# print(coso)
#
# coso2 = pn.getAtomDict(['Ar3', 'Ar4', 'He1'], only_coll=True)
#
# print(coso2)

