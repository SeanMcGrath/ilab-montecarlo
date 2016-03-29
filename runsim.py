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

ref_data = np.array([0, 0, 4, 7, 13, 15, 4])
sum_residuals = []

for muon_mass in range(60, 140, 5):
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
	hist = np.array([(sum(ref_data)/num_decays)*el for el in hist])
	residuals = np.square(ref_data - hist)
	print(residuals.sum())
	sum_residuals.append((muon_mass, residuals.sum()))
	# plt.hist(sparks, bins=[el-0.5 for el in range(1, 9)])
	# plt.xlim(-0.5, 7.5)
	# plt.title('spark counts for muon mass = ' + str(muon_mass))
	# plt.show()

print(min(sum_residuals, key=lambda x: x[1]))