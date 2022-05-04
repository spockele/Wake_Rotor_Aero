from logging import root
import numpy as np

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
        if not bLocalCylindrical:
            # Input is in carthesian
            self.xloc = localPos[0]
            self.yloc = localPos[1]
            self.zloc = localPos[2]
        
        if bLocalCylindrical:
            self.xloc, self.yloc, self.zloc = CylToCarth(localPos)

        # Setting the global coordinates
        self.xglob = self.xloc + referencePos[0]
        self.yglob = self.yloc + referencePos[1]
        self.zglob = self.zloc + referencePos[2]

        # Now update the cylindrical positions
        self.rloc, self.thetaloc, _ = CarthToCyl(self.xloc, self.yloc, self.zloc)

    def Length(self):
        return np.sqrt(self.xloc*self.xloc + self.yloc*self.yloc + self.zloc*self.zloc)

class ControlPoint:
    '''Position of control point in reference frame (which has root in rootPos)'''
    def __init__(self, startPos, endPos = None):
        self.startpos = startPos

        self.circulation = None
        self.endPoint = None

    def set_circulation(self, magnitude):
        self.circulation = magnitude



class Leg:
    def __init__(self, reset=False):
        if not reset:
            self.control_points = []

        self.circulation = None

    def SetControlPoints(self, arrayOfCylindricalPositions, rootPos):
        '''Takes an argument like arrayOfCylindricalPositions = [(2,pi,0), (2.1,pi,1),...] and a root reference position in carthesian.'''

        for cylPos in arrayOfCylindricalPositions:
            point = Vec(cylPos, rootPos, bLocalCylindrical=True)
            # Set the previous point to end in this new control point
            self.control_points[-1].endPoint = point
            self.control_points.append(ControlPoint(point))
        

            

    def reset(self):
        [cp.reset() for cp in self.control_points]
        self.__init__(reset=True)


class HorseShoe:
    def __init__(self, reset=False):
        if not reset:
            self.leg1 = Leg()
            self.leg2 = Leg()

        self.circulation = None

    def reset(self):
        self.leg1.reset()
        self.leg2.reset()
        self.__init__(reset=True)


class Wake:
    def __init__(self, reset=False):
        if not reset:
            self.horse_shoes = []

    def reset(self):
        [hs.reset() for hs in self.horse_shoes]
        self.__init__(reset=True)

if __name__=="__main__":
    print("This is a lifting line library, pls dont run this")