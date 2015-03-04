import os
from math import *
from mks import *
import numpy as np

def elements_radians2degrees(elements_in):
    elements = elements_in.copy()
    elements[:,2:] = elements[:,2:]*180/pi
    return(elements)


def elements_to_cartesian(elements):
    """Convert orbial elements to cartesian coordinates.

    Args:
        elements: (nx6) array where each row gives the orbital elements for an orbiting body.
    Returns:
        (nx6) array of cartesian coordinates

    """
    coordinates = np.array([])
    for ii in range(len(elements)):
        a = elements[ii,0]
        e = elements[ii,1]
        i = elements[ii,2]
        O = elements[ii,3]
        o = elements[ii,4]
        M = elements[ii,5]
        #1. calculate the eccentric anomaly
        E = M
        if e != 0 and M != 0:
            Enext = e*sin(E) + M
            while (abs(Enext-E)/Enext > 1e-12):
                print "            ", (abs(Enext-E)/Enext)
                E = Enext
                Enext = e*sin(E) + M
            E = Enext
        #2. Calculate coords in orbital plane reference
        n = a**(-3./2) # GM=1
        q = np.array([a*(cos(E)-e), a*sqrt(1-e*e)*sin(E), 0])
        qdot = np.array([-sin(E), sqrt(1-e*e)*cos(E), 0]) * n*a/(1-e*cos(E))

        #3. Now rotate into reference frame.
        r3O = np.array([[cos(-O), sin(-O), 0], [-sin(-O), cos(-O), 0], [0, 0, 1.]])
        r1i = np.array([[1., 0, 0], [0, cos(-i), sin(-i)], [0, -sin(-i), cos(-i)]])
        r3o = np.array([[cos(-o), sin(-o), 0], [-sin(-o), cos(-o), 0], [0, 0, 1.]])
        r = np.dot(r3O, np.dot(r1i, np.dot(r3o, q)))
        rdot = np.dot(r3O, np.dot(r1i, np.dot(r3o, qdot)))
        coordinates = np.hstack((coordinates, r, rdot))
    coordinates.shape = (len(coordinates)/6,6)
    return(coordinates)

def cartesian_to_elements(coords):
    """Convert cartesian coordinates to orbital elements.

    Args:
        coords: (nx6) array of cartesian coordinates
    Returns:
        (nx6) array where each row gives the orbital elements for an orbiting body.

    """

    elements = np.array([])
    for ii in range(len(coords)):
        r = coords[ii,:3]
        rdot = coords[ii,3:]
        h = np.cross(r, rdot)
        O = atan2(h[0], -h[1])
        i = atan2(sqrt(h[0]*h[0]+h[1]*h[1]), h[2])
        R1i = np.array([[1., 0, 0], [0, cos(i), sin(i)], [0, -sin(i), cos(i)]])
        R3O = np.array([[cos(O), sin(O), 0], [-sin(O), cos(O), 0], [0, 0, 1.]])
        p = np.dot(R1i, np.dot(R3O,r))
        u = atan2(p[1], p[0])
        hmag = sqrt(np.dot(h,h))
        rmag = sqrt(np.dot(r,r))
        vmag = sqrt(np.dot(rdot, rdot))
        a = rmag/(2-rmag*vmag*vmag)
        e = sqrt(abs(1-hmag*hmag/a))
        if e < 1e-6:
            e=0
        radial_vel = np.dot(r, rdot)/rmag
        if e != 0:
            try:
                E = acos((a-rmag)/a/e)
            except ValueError:
                if (a-rmag)/a/e > 1:
                    E = 0
                elif (a-rmag)/a/e < -1:
                    E = pi
            if abs(E - asin(rmag*radial_vel/e/sqrt(a)))>pi/2:
                E = -E
            nu = atan2(sqrt(1-e*e)*sin(E), cos(E)-e)
        else:
            E = u
            nu = u
        o = u - nu
        M = E-e*sin(E)
        elements = np.hstack((elements, np.array([a, e, i, O, o, M])))
    elements.shape = (len(elements)/6, 6)
    return(elements)



def transit_observable(orbits, reference, r_star=1):
    """ Given the orbital elements will the planet transit from
    point of view of reference direction.

    Args:
        orbits: Orbital elements of the bodies
            a: Semi-major axis [AU]
            e: eccentricity
            i: inclination [radians]
            Omega: Longitude of the ascending node [radians]
            omega: Argument of periapse [radians]
            m: mean anomaly (not used) [radians]
        r_star: Stellar radius [R_sun]

    Return:
        true or false whether the planet will transit.

    """
    orbits[:,5] = 0
    cartesian = elements_to_cartesian(orbits)
    transits = []


    for ii in range(len(orbits)):
        a = orbits[ii,0]
        e = orbits[ii,1]
        # Transform the reference direction into the coordinate system of the orbit
        x_axis = cartesian[ii,:3].copy()
        z_axis = np.cross(cartesian[ii,:3], cartesian[ii,3:])
        y_axis = np.cross(z_axis, x_axis)
        x_axis = x_axis/np.linalg.norm(x_axis)
        y_axis = y_axis/np.linalg.norm(y_axis)
        z_axis = z_axis/np.linalg.norm(z_axis)
        rot = np.array([x_axis, y_axis, z_axis])
        reference = np.dot(rot, reference)

        # Calculate the impact parameter
        nu = atan2(reference[1], reference[0])
        orbital_distance_at_nu = a*(1-e*e)/(1+e*cos(nu))
        height_above_orbital_plane = reference[2]/sqrt(np.dot(reference[:2], reference[:2]))*orbital_distance_at_nu
        reference_angle = atan(reference[2]/np.linalg.norm(reference[:2]))
        projected_height = height_above_orbital_plane*cos(reference_angle)
        transits.append(abs(projected_height*AU/RSUN) < r_star)
    return(transits)


class MercuryInputWriter:
    def __init__(self):
        self.smallContent = ")0+_06 Small-body initial data (WARNING: Do not"\
                            " delete this line!!)\n)------------------------"\
                            "-------------------------------\nstyle = Astero"\
                            "id\n"
        self.bigContent = ")0+_06 Small-body initial data (WARNING: Do not delete this line!!)\n)-------------------------------------------------------\nstyle = Asteroid\nepoch=0\n"
        self.directory = "./"

    def AddBigBody(self, body):
        self.bigContent += body.InputFormattedString()

    def AddSmallBody(self, body):
        self.smallContent += body.InputFormattedString()

    def WriteFiles(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)
        smallFile = open(self.directory+"small.in", "w")
        bigFile = open(self.directory+"big.in", "w")
        smallFile.write(self.smallContent)
        bigFile.write(self.bigContent)
        smallFile.close()
        bigFile.close()

    def ClearBodies(self):
        self.smallContent = ")0+_06 Small-body initial data (WARNING: Do not"\
                            " delete this line!!)\n)------------------------"\
                            "-------------------------------\nstyle = Astero"\
                            "id\n"
        self.bigContent = ")0+_06 Small-body initial data (WARNING: Do not delete this line!!)\n)-------------------------------------------------------\nstyle = Asteroid\nepoch=0\n"


class MercuryBody:
    def __init__(self, name, a, e, i, g, n, m, mass=None, r=None, density=None,
             a1=None, a2=None, a3=None, ep=None):
        """Inputs:
        a = semi-major axis (in AU)
        e = eccentricity
        i = inclination (degrees)
        g = argument of pericentre (degrees)
        n = longitude of the ascending node (degrees)
        m = mean anomaly (degrees)

        ep = epoch of osculation

        mass = X    where X is a real number, to indicate the body's mass in
        Solar masses. If you don't specify a value the mass is
        assumed to be 0.

        r = X    where X is a real number, to indicate the maximum
        distance from the body (in Hill radii) that constitutes
        a close encounter. If you don't include this the default
        is r=1

        density = X    where X is a real number, to indicate the density of the
        body in g/cm^3. If you don't include this the default is d=1

        a1 = X   where X is a real number, to indicate the A1 non-gravitational
        force parameter for this body. Realistically this should be
        zero for Big bodies (the default is 0).

        a2 = X   where X is a real number, to indicate the A2 non-gravitational
        force parameter for this body (the default is 0).

        a3 = X   where X is a real number, to indicate the A1 non-gravitational
        force parameter for this body (the default is 0).
        """
        self.name=name
        self.orbitalElements=[a,e,i,g,n,m]
        self.bodyProperties={'m':mass, 'r':r, 'd':density, 'a1':a1, 'a2':a2,
                             'a3':a3, 'ep':ep}
        for key in self.bodyProperties.keys():
            if self.bodyProperties[key] == None:
                del self.bodyProperties[key]

    def InputFormattedString(self):
        string = self.name
        for key in self.bodyProperties.keys():
            string += "  "+key+"="+str(self.bodyProperties[key])
        string += "\n"
        for element in self.orbitalElements:
            string += "  "+str(element)
        string += "  0  0  0\n"
        return string