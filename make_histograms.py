#!/bin/env/ python

import sys
import math
import random
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt

from resources import *

sns.set(context='poster', font_scale=1.5, style='white')

try:
	num_decays = int(sys.argv[1])
except:
	num_decays = 10000

ref_x = [1, 2, 3, 4, 5, 6, 7]
ref_data = np.array([0, 0, 4, 7, 13, 15, 4])

masses = [90, 100]

f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)

for muon_mass, ax in zip(masses, (ax1, ax2)):
	# generate random energies
	energy_pdf = []
	energies = []
	while len(energies) < num_decays:
		energy = random.uniform(0, muon_mass/2)
		energy_prob  = random.uniform(0, fermi(muon_mass, muon_mass/2))
		if energy_prob < fermi(muon_mass, energy):
			energy_pdf.append(energy_prob)
			energies.append(energy)

	# generate random r values
	possible_r = np.arange(0, R/2, 0.001)
	r_rand = GeneralRandom(possible_r, possible_r)
	rs = [(R/2)*el for el in r_rand.random(1000000) if el != 0]

	# run simulation
	electrons = []
	for energy, r in zip(energies, rs):
		electrons.append(Electron(muon_mass, energy, r))

	# get spark counts
	sparks = [e.sparks for e in electrons]

	hist = np.histogram(sparks, bins=[el-0.5 for el in range(1, 9)])[0]
	scaled_hist = np.array([(sum(ref_data)/num_decays)*el for el in hist])

	# plt.hist(sparks, bins=[el-0.5 for el in range(1, 9)])
	ax.bar(ref_x, scaled_hist, label='Theoretical', color='0.75')
	ax.errorbar([el + .4 for el in ref_x], ref_data, yerr=np.sqrt(ref_data), fmt='ko', label='Experimental')

ax1.set_title('Spark Counts for Different Muon Masses')
ax1.text(7, 15, '90 MeV')
ax2.text(7, 15, '100 MeV')
ax1.legend(loc=2)
plt.show()
