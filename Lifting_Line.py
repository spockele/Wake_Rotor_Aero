import numpy as np


class ControlPoint:
    def __init__(self, r, theta, z, x0, y0):
        self.r = r
        self.theta = theta
        self.z = z

        self.x = r * np.cos(theta) + x0
        self.y = r * np.sin(theta) + y0

        self.circulation = None
        self.orientation = None

    def set_circulation(self, magnitude, orientation):
        self.circulation = magnitude
        self.orientation = orientation


class Leg:
    def __init__(self):
        self.control_points = []

        self.circulation = None


class HorseShoe:
    def __init__(self):
        self.leg1 = Leg()
        self.leg2 = Leg()

        self.circulation = None


class Wake:
    def __init__(self):
        self.horse_shoes = []
