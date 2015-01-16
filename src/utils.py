import os

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