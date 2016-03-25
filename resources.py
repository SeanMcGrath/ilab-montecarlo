import random
import math
import functools
import numpy as np

# constants
CHI_0 = 8.9 #cm
E_CRIT = 45 #MeV
R = 7.62 #cm
PLATE_THICKNESS = 0.9525 #cm
GAP_THICKNESS = 0.625 #cm

# Generate random val of m_mu_c^2 in MeV
def random_m_mu():
    return random.randint(60,140)

def fermi(muon_mass, energy):
    return (muon_mass * energy)**2 * (3 - 4*energy/muon_mass)


class GeneralRandom:
  """This class enables us to generate random numbers with an arbitrary 
  distribution."""
  
  def __init__(self, x = np.arange(-1.0, 1.0, .01), p = None, Nrl = 1000):
    """Initialize the lookup table (with default values if necessary)
    Inputs:
    x = random number values
    p = probability density profile at that point
    Nrl = number of reverse look up values between 0 and 1"""  
    if p is None:
      p = np.exp(-10*x**2.0)
    self.set_pdf(x, p, Nrl)
    
  def set_pdf(self, x, p, Nrl = 1000):
    """Generate the lookup tables. 
    x is the value of the random variate
    pdf is its probability density
    cdf is the cumulative pdf
    inversecdf is the inverse look up table
    
    """
    
    self.x = x
    self.pdf = p/p.sum() #normalize it
    self.cdf = self.pdf.cumsum()
    self.inversecdfbins = Nrl
    self.Nrl = Nrl
    y = np.arange(Nrl)/float(Nrl)
    delta = 1.0/Nrl
    self.inversecdf = np.zeros(Nrl)    
    self.inversecdf[0] = self.x[0]
    cdf_idx = 0
    for n in range(1,self.inversecdfbins):
      while self.cdf[cdf_idx] < y[n] and cdf_idx < Nrl:
        cdf_idx += 1
      self.inversecdf[n] = self.x[cdf_idx-1] + (self.x[cdf_idx] - self.x[cdf_idx-1]) * (y[n] - self.cdf[cdf_idx-1])/(self.cdf[cdf_idx] - self.cdf[cdf_idx-1]) 
      if cdf_idx >= Nrl:
        break
    self.delta_inversecdf = np.concatenate((np.diff(self.inversecdf), [0]))
              
  def random(self, N = 1000):
    """Give us N random numbers with the requested distribution"""

    idx_f = np.random.uniform(size = N, high = self.Nrl-1)
    idx = np.array([idx_f],'i')
    y = (self.inversecdf[idx] + (idx_f - idx)*self.delta_inversecdf[idx])[0]

    return y
  

class Electron:
    
    Z0_LIMITS = (0, PLATE_THICKNESS) #cm
    THETA_LIMITS = (0, math.pi/6) #rad
    PHI_LIMITS = (0, 2*math.pi) #rad

    def __init__(self, muon_mass, energy, r):
        self.z0     = random.uniform(*self.Z0_LIMITS)
        self.theta  = random.uniform(*self.THETA_LIMITS)
        self.phi    = random.uniform(*self.PHI_LIMITS)
        self.energy = energy
        self.r      = r

    def __repr__(self):
        return "track: {} escape: {} sparks: {}".format(
            self.track_length, self.escape_length, self.sparks)
    
    @property
    def track_length(self):
        return CHI_0 * math.log(1 + self.energy/E_CRIT)

    @property
    def escape_length(self):
        numerator = -self.r * math.cos(self.phi) + \
                math.sqrt(self.r**2 * math.cos(self.phi)**2 + (R**2 - self.r**2))
        denominator = math.sin(self.theta)
        frac = PLATE_THICKNESS / (PLATE_THICKNESS + GAP_THICKNESS)

        return (numerator / denominator) * frac

    @property
    def sparks(self):
        track = self.track_length
        escape = self.escape_length
        if track > escape:
            track = escape
        return 1 + math.floor((track * math.cos(self.theta) - self.z0)/PLATE_THICKNESS)

    def weighted_sparks(self, muon_mass):
        return self.r * (muon_mass * self.energy)**2 * (3 - 4*self.energy/muon_mass) * self.sparks
    