import numpy as np

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