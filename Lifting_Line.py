from logging import root
from multiprocessing.sharedctypes import Value
import numpy as np
from BEM_code import DU95W150


"""
ALL ANGLES IN RADIANS!!!
"""

# Cylindrical position vector: (r, theta=radians, z)
# Carthesian position vector: (x, y, z)

def CylToCarth(r,theta,z):
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return (x,y,z)
    
def CarthToCyl(x,y,z):
    r = np.sqrt(x*x + y*y)
    theta = np.arctan2(y,x)
    return (r, theta, z)

class Vec:
    '''Reference pos is always in carthesian coordinates.'''
    def __init__(self,localPos, referencePos, bLocalCylindrical=False):
        self.refPos = referencePos
        if not bLocalCylindrical:
            # Input is in carthesian
            self.xloc = localPos[0]
            self.yloc = localPos[1]
            self.zloc = localPos[2]
        
        if bLocalCylindrical:
            self.xloc, self.yloc, self.zloc = CylToCarth(*localPos)
            self.rloc, self.thetaloc, _ = localPos

        # Setting the global coordinates
        self.xglob = self.xloc + referencePos[0]
        self.yglob = self.yloc + referencePos[1]
        self.zglob = self.zloc + referencePos[2]

    def Length(self):
        return np.sqrt(self.xloc*self.xloc + self.yloc*self.yloc + self.zloc*self.zloc)

    # Operator overloading

    def __repr__(self):
        return "carth=({:.3f}, {:.3f}, {:.3f}), cyl=({:.3f}, {:.3f}, {:.3f}), ref=({:.3f}, {:.3f}, {:.3f})".format(self.xloc, self.yloc, self.zloc, self.rloc, self.thetaloc, self.zloc, self.refPos[0], self.refPos[1], self.refPos[2])

    def __add__(self, o):
        # Note that order of operation matters: The reference position of the FIRST entry is kept!
        x = self.xloc + (o.xglob - self.refPos[0])
        y = self.yloc + (o.yglob - self.refPos[1])
        z = self.zloc + (o.zglob - self.refPos[2])
        return Vec((x,y,z), self.refPos)

    def __mul__(self, scale):
        x = self.xloc * scale
        y = self.yloc * scale
        z = self.zloc * scale
        return Vec((x,y,z), self.refPos)

    def __sub__(self, o):
        return self + -1. * o

    def __div__(self, scale):
        return self * (1/scale)



class Filament:
    ''' A Vortex Filament from startpos to endpos with strength circulation. '''
    def __init__(self, startPos, endPos):
        self.startPos = startPos
        self.endPos = endPos
        self.centre = (endPos - startPos)/2

        self.circulation = None

    def set_circulation(self, magnitude):
        self.circulation = magnitude



class Leg:
    def __init__(self, reset=False):
        if not reset:
            self.control_points = []

        self.circulation = None

    def CreateControlPoints(self, arrayOfCylindricalPositions, rootPos):
        '''Takes an argument like arrayOfCylindricalPositions = [(2,pi,0), (2.1,pi,1),...] and a root reference position in carthesian.'''
        
        # Requires a starting point!
        prevPoint = Vec(arrayOfCylindricalPositions[0])
        arrayOfCylindricalPositions.pop(0)
        for cylPos in arrayOfCylindricalPositions:
            point = Vec(cylPos, rootPos, bLocalCylindrical=True)
            self.control_points.append(Filament(point, prevPoint))  
            prevPoint = point     

    def reset(self):
        [cp.reset() for cp in self.control_points]
        self.__init__(reset=True)


class HorseShoe:
    def __init__(self, airfoil, chord, delta_r, twist, pos, reset=False):
        if not reset:
            self.leg1 = Leg()
            self.leg2 = Leg()

            self.airfoil = airfoil
            self.chord = chord
            self.delta_r = delta_r
            self.twist = twist
            self.pos = pos

        self.circulation = None
        self.induced_velocity = None

    def set_circulation(self, w, aoa):
        """
        Determine the circulation based on the inflow conditions and the airfoil lift curve
        :param w: inflow velocity
        :param aoa: angle of attack
        :return: None
        """
        self.circulation = .5 * w * self.delta_r * self.chord * self.airfoil.cl(np.degrees(aoa))

    def set_velocity_triangle(self):
        return

    def reset(self):
        self.leg1.reset()
        self.leg2.reset()
        self.__init__(reset=True)


class Wake:
    """
    INPUTS: N is the number blade elements
    """
    def __init__(self,N , reset=False):
        if not reset:
            self.N = N
            self.horse_shoes = []

    def reset(self):
        [hs.reset() for hs in self.horse_shoes]
        self.__init__(reset=True)

if __name__=="__main__":
    print("This is a lifting line library, pls dont run this")