#!/bin/env/ python

import sys
import math
import random
import seaborn as sns
import numpy as np
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

energies = [el for el in fermi_rand.random(1000000) if el != 0]

plt.hist(energies)
plt.show()
plt.hist(rs)
plt.show()