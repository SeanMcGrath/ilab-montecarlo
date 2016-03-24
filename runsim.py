#!/bin/env/ python

import sys
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt

from resources import *

try:
	num_decays = int(sys.argv[1])
except:
	num_decays = 100000

for mass in [100]:

	electrons = [Electron(mass) for i in range(num_decays)]

	sparks = [e.sparks for e in electrons]

	plt.hist(sparks)
	plt.title('spark counts for muon mass = ' + str(mass))
	plt.show()