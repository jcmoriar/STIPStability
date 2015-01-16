import numpy as np
from math import *
from utils import *
import matplotlib.pyplot as plt
from mks import *

PMIN = 1e0
PMAX = 200e0
MINMASS = 5e0
MAXMASS = 10e0
MSTAR = 0.4


#constants
DAYSPERYEAR = 365.242


class PlanetarySystem(object):
    """ Orbital data for a planetary system.

    Currently, the system is constructed with initial eccentricities
    of 0, inclinations following a rayleigh distribution and orbital
    periods distributed uniformly in log space.

    """

    def __init__(self, num_bodies, sigma_inc):
        self.star_mass = MSTAR
        self.periods = np.random.rand(num_bodies)*(log10(PMAX)-log10(PMIN))+log10(PMIN)
        self.periods = (np.power(10., self.periods))
        self.periods.sort()
        self.masses = np.random.rand(num_bodies)*(MAXMASS - MINMASS) + MINMASS
        self.inclinations = np.random.rayleigh(sigma_inc, num_bodies)
        self.sma = np.power(np.power(self.periods/DAYSPERYEAR, 2)*self.star_mass, 1./3.)
        self.eccentricities = np.zeros(len(self.periods))

    def write_mercury_input(self, dir="./"):
        writer = MercuryInputWriter()
        writer.directory = dir
        for i in range(len(self.periods)):
            planet = MercuryBody("P"+str(i), self.sma[i], self.eccentricities[i],
                                 self.inclinations[i], np.random.rand()*360.,
                                 np.random.rand()*360., np.random.rand()*360.,
                                 mass=self.masses[i]*MEARTH/MSUN, density=3e0)
            writer.AddBigBody(planet)
        writer.WriteFiles()


    def plot(self):
        fig = plt.figure()
        axis = fig.add_subplot(111)
        axis.semilogx(self.periods, self.mass)
