# Grotrian diagram of [O III]

# Imports
import pyneb as pn
from matplotlib import pyplot as plt
# # Define the atom

n2 = pn.Atom('N', 2)

tem,den = 10000, 1000
print(n2.getEmissivity(tem, den, wave=6583)/n2.getEmissivity(tem, den, wave=6527))

for atom, ionization in zip(('O', 'S', 'N'), (3, 2, 2)):
    ion = pn.Atom(atom, ionization)

    # Draw the diagram
    ion.plotGrotrian(detailed=True, thresh_int=1e-5)
    plt.show()
