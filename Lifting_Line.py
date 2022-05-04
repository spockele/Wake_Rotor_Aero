import numpy as np


class ControlPoint:
    def __init__(self, r, theta, z, x0, y0):
        self.r = r
        self.theta = theta
        self.z = z

        self.x0 = x0
        self.y0 = y0

        self.x = r * np.cos(theta) + x0
        self.y = r * np.sin(theta) + y0

        self.circulation = None
        self.orientation = None

    def set_circulation(self, magnitude, orientation):
        self.circulation = magnitude
        self.orientation = orientation

    def reset(self):
        self.__init__(self.r, self.theta, self.z, self.x0, self.y0)


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