import random
import math
import functools

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

# define any function here!
def random_from_distribution(distribution_func, xmin, xmax):
    # f(x) = 1.0 : for uniform probability distribution

    # f(x) = x : for triangular probability distribution
    # (math.sqrt(random.random()) would also produce triangular p.d. though.)

    # f(x) = math.exp(-x*x/2.0)/math.sqrt(2.0*math.pi) : for std normal p.d.
    # (taking average of (last) 2,3,... random.random() values would also
    # produce normal probability distributions though.)

    to_return = []

    # find ymin-ymax
    numSteps = 1000000 # bigger the better but slower!
    ymin = distribution_func(xmin)
    ymax = ymin
    for i in range(numSteps):
        x = xmin + (xmax - xmin) * float(i) / numSteps
        y = distribution_func(x)
        if y < ymin: ymin = y
        if y > ymax: ymax = y

    while y <= distribution_func(x):
        # generate a random number between 0 to 1
        xr = random.random()
        yr = random.random()
        x = xmin + (xmax - xmin) * xr
        y = ymin + (ymax - ymin) * yr
    
    return y


class Electron:
    
    Z0_LIMITS = (0, PLATE_THICKNESS) #cm
    R_LIMITS = (0, R/2) #cm
    THETA_LIMITS = (0, math.pi/6) #rad
    PHI_LIMITS = (0, 2*math.pi) #rad

    def __init__(self, muon_mass):
        self.z0     = random.uniform(*self.Z0_LIMITS)
        self.r      = random.uniform(*self.R_LIMITS)
        self.theta  = random.uniform(*self.THETA_LIMITS)
        self.phi    = random.uniform(*self.PHI_LIMITS)
        self.energy = random_from_distribution(functools.partial(fermi, muon_mass), 0, muon_mass/2)

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
    