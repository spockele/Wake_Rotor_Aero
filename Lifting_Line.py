import numpy as np
from BEM_code import DU95W150

def read_from_file(path):
    f = open(path)
    lines = f.readlines()
    out_list = [[float(num) for num in line.strip('\n').split(',')] for line in lines]
    return np.array(out_list)

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
    def __init__(self, airfoil, chord, r_inner, r_outer, relative_pitch, reset=False):
        self.delta_r = r_outer-r_inner #length of lifting line element or blade element
        if not reset:
            self.leg1 = Leg()
            self.leg2 = Leg()

            self.airfoil = airfoil
            self.chord = chord
            self.delta_r = delta_r

        self.circulation = None

    def set_circulation(self, w, aoa):
        """
        Determine the circulation based on the inflow conditions and the airfoil lift curve
        :param w: inflow velocity
        :param aoa: angle of attack
        :return: None
        """
        self.circulation = .5 * w * self.delta_r * self.chord * self.airfoil.cl(np.degrees(aoa))

    def reset(self):
        self.leg1.reset()
        self.leg2.reset()
        self.__init__(reset=True)


class Turbine:
    """
    Turbine parameters, air density, U_inf
    """
    def __init__(self , reset=False):
        data = read_from_file('DU95W150.csv')
        self.alpha_lst = data[:, 0]
        self.cl_lst = data[:, 1] #; self.cd_lst = data[:, 2]; self.cm_lst = data[:, 3]

        self.b = 3 # Number of blades
        self.n_elements = 1 # Divide the blade up in n_elements
        self.rho = 1.225 # [Kg/m3]
        self.u_inf = 10 # [m/s] U_infinity = free stream velocity
        self.radius = 50 # Total radius
        self.blade_pitch = 0
        r_start = 0.2*self.radius

        self.bladeElement = list()

        for i in range(self.n_elements):
            r_inner = r_start + (self.radius - r_start) / self.n_elements * i
            r_outer = r_start + (self.radius - r_start) / self.n_elements * (i+1)
            r = np.mean([r_inner, r_outer])
            # Sorry for hardcoding the equations below- taken from the assignment description :)
            twist = 14 * (1 - r / self.radius)
            chord = (3 * (1 - r / self.radius) + 1)

            # BladeElement takes in argument relative_pitch, I assume that this means total? So offset with the blade pitch
            relative_pitch = self.blade_pitch + twist

            self.bladeElement.append(chord, r_inner, r_outer, relative_pitch)

        if not reset:
            self.horse_shoes = []

    def cl(self, alpha):
        return np.interp(alpha, self.alpha_lst, self.cl_lst)

    def reset(self):
        [hs.reset() for hs in self.horse_shoes]
        self.__init__(reset=True)




if __name__=="__main__":
    print("This is a lifting line library, pls dont run this")
    Turbine()