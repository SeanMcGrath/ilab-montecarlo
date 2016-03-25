#!/bin/env/ python

import sys
import seaborn as sns
from matplotlib import pyplot as plt

from resources import *

try:
	num_decays = int(sys.argv[1])
except:
	num_decays = 10000

muon_mass = 100

possible_energies = np.arange(0, muon_mass/2, .01)
energy_pdf = [fermi(muon_mass, energy) for energy in possible_energies]

fermi_rand = GeneralRandom(possible_energies, np.array(energy_pdf))

possible_r = np.arange(0, R/2, 0.001)

r_rand = GeneralRandom(possible_r, possible_r)

rs = [(R/2)*el for el in r_rand.random(1000000) if el != 0]

energies = [5*el for el in fermi_rand.random(1000000) if el != 0]
electrons = []
for i in range(num_decays):
	electrons.append(Electron(muon_mass, energies[i], rs[i]))

sparks = [e.sparks for e in electrons]

plt.hist(sparks, bins=[el-0.5 for el in range(9)])
plt.xlim(-0.5, 7.5)
plt.title('spark counts for muon mass = ' + str(muon_mass))
plt.show()