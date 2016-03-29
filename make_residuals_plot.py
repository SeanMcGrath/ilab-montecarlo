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

ref_data = np.array([0, 0, 4, 7, 13, 15, 4])
sum_residuals = []
masses = range(60, 145, 5)

for muon_mass in masses:
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
	sum_residuals.append(residuals.sum())

plt.plot(masses, sum_residuals, 'ko')
plt.plot(masses, sum_residuals, 'k--')
plt.ylabel('Sum of Squared Residuals')
plt.xlabel('Muon Mass (MeV)')
plt.xlim(50, 150)
plt.annotate('Minimum', xy=(masses[8], sum_residuals[8] + 5), xytext=(masses[8] - 5, sum_residuals[8] + 45),\
	arrowprops=dict(facecolor='black', shrink=0.05))
plt.title('Effect of Model Muon Mass on Model Residuals')
plt.show()