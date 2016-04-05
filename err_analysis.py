#!/bin/env/ python

import sys
import math
import random
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from lmfit import Parameters, conf_interval, minimize, fit_report

from resources import *

try:
	num_decays = int(sys.argv[1])
except:
	num_decays = 10000


def add1(x):
	if x == 0:
		return 1
	else:
		return x

ref_data = np.array([0, 0, 4, 7, 13, 15, 4])
ref_err = np.array(list(map(add1, np.sqrt(ref_data))))
possible_r = np.arange(0, R/2, 0.001)

def residuals(parameters):
	# generate random energies
	muon_mass = parameters['mass'].value
	energy_pdf = []
	energies = []
	while len(energies) < num_decays:
		energy = random.uniform(0, muon_mass/2)
		energy_prob  = random.uniform(0, fermi(muon_mass, muon_mass/2))
		if energy_prob < fermi(muon_mass, energy):
			energy_pdf.append(energy_prob)
			energies.append(energy)

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

	residuals = ref_data - hist
	weighted = np.sqrt(residuals**2 / ref_err**2)

	return weighted

params = Parameters()
params.add('mass')
chis = []
masses = range(80,120)

for m in masses:
	params['mass'].set(m)
	mi = minimize(residuals, params)
	chis.append(mi.redchi)

plt.plot(masses, chis, 'ko')
plt.plot(masses, chis, 'k--')
plt.show()