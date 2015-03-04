import utils as u
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *


def mutually_observable_planets(orbits, num_trials=10000):
    """ Calculate the number of observable transits for a system
    from num_trials different viewing angles.

    Angles units should be in degrees.
    """

    orbits[:,2:] = orbits[:,2:]*pi/180
    num_observable=[]
    for i in range(num_trials):
        # Create a randomly distributed reference direction
        sign = np.random.rand()*2-1
        sign = sign/abs(sign)
        theta = acos(2*np.random.rand()-1)*sign
        phi = 2*pi*np.random.rand()
        reference_diretion = np.array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])
        # See how many transits are observable from that direction
        transits = u.transit_observable(orbits, reference_diretion)
        num_observable.append(transits.count(True))
    # print num_observable.count(1)/float((num_observable.count(0)+num_observable.count(1)))
    return(num_observable)

def planets_observed_hist(num_observable_in):
    num_observable = []
    for num in num_observable_in:
        if num !=0:
            num_observable.append(num)
    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.hist(num_observable, bins=[1,2,3,4,5,6,7,8,9,10], align="left")
    axis.set_xlabel("Number of planets observed")




