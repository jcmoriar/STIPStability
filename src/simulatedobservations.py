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
        transits = u.transit_observable(orbits)
        num_observable.append(transits.count(True))
    print num_observable.count(1)/float((num_observable.count(0)+num_observable.count(1)))
    return(num_observable)

def planets_observed_hist(num_observable):
    fig = plt.figure()
    axis = fig.add_subplot(111)
    axis.hist(num_observable)
    axis.set_xlabel("Number of planets observed")

def plot_angles(x, y, z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x,y,z)
    ax.set_xlabel("x")
    ax.set_ylabel("y")

def test_rot(orbit, axis, theta):
    orbit = np.array(orbit)
    orbit[2:] = orbit[2:]*pi/180
    orbit.shape=(1,6)
    axis = np.array(axis)
    theta = theta*pi/180
    new_orbit = rotate_system(axis, theta, orbit)
    new_orbit[0,2:] = new_orbit[0,2:]/pi*180
    print new_orbit.tolist()


