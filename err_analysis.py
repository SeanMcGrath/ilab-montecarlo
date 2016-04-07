#!/bin/env/ python

import sys
import math
import random
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import minimize

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

def chisq(ydata,ymod,sd=None):  

      # Chi-square statistic (Bevington, eq. 6.9)  
      if sd is None:  
           chisq=np.sum((ydata-ymod)**2)  
      else:  
           chisq=np.sum( ((ydata-ymod)/sd)**2 )  
        
      return chisq  

ref_x = np.arange(1,8)
ref_data = np.array([0, 0, 4, 7, 13, 15, 4])
ref_err = np.array(list(map(add1, np.sqrt(ref_data))))
possible_r = np.arange(0, R/2, 0.001)

def monte(muon_mass):
	# generate random energies
	print(muon_mass)
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
	hist = hist * sum(ref_data)/num_decays
	return chisq(ref_data, hist, ref_err)

results = []
for i in range(2):
	res = minimize(monte, x0=(95), method='powell', options={'disp': True})
	print(res.x)
	results.append(res.x)

print(np.mean(results))
print(np.std(results))
