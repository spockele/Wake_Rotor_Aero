import numpy as np
from BEM_code import DU95W150


"""
ALL ANGLES IN RADIANS!!!
"""

# Cylindrical position vector: (r, theta, z)
# Carthesian position vector: (x, y, z)

class ControlPoint:
    '''Position of control point in reference frame (which has root in rootPos)'''
    def __init__(self, cylPos, rootPos):
        self.cylPos = cylPos

        self.origin = rootPos

        x = cylPos(0) * np.cos(cylPos(1)) + rootPos(0)
        y = cylPos(0) * np.sin(cylPos(1)) + rootPos(1)
        self.pos = (x,y)

        self.circulation = None
        self.orientation = None

    def set_circulation(self, magnitude, orientation):
        self.circulation = magnitude
        self.orientation = orientation



class Leg:
    def __init__(self, reset=False):
        if not reset:
            self.control_points = []

        self.circulation = None

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