from logging import root
from multiprocessing.sharedctypes import Value
import numpy as np
from BEM_code import DU95W150
import matplotlib.pyplot as plt


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
    def __init__(self,localPos, referencePos=(0.,0.,0.), bLocalCylindrical=False):

        # I make this error too much so this is an attempt at salvation
        if bLocalCylindrical < 0 or bLocalCylindrical > 1:
            raise ValueError("b-baka, are you sure you passed a tuple into the constructor for Vec like Vec((1,2,3)), and not Vec(1,2,3)?")

        self.refPos = referencePos
        if not bLocalCylindrical:
            # Input is in carthesian
            self.xloc, self.yloc, self.zloc = localPos
            self.rloc, self.thetaloc, _ = CarthToCyl(*localPos)
        else:
            self.xloc, self.yloc, self.zloc = CylToCarth(*localPos)
            self.rloc, self.thetaloc, _ = localPos

        # Setting the global coordinates
        self.xglob = self.xloc + referencePos[0]
        self.yglob = self.yloc + referencePos[1]
        self.zglob = self.zloc + referencePos[2]

    def Length(self):
        return np.sqrt(self.xloc*self.xloc + self.yloc*self.yloc + self.zloc*self.zloc)

    def output_xyz(self, glob=False):
        x = self.xglob if glob else self.xloc
        y = self.yglob if glob else self.yloc
        z = self.zglob if glob else self.zloc

        return x, y, z

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

    def __rmul__(self, scale):
        return self * scale

    def __sub__(self, o):
        return self + -1. * o

    def __truediv__(self, scale):
        return self * (1/scale)


class Filament:
    ''' A Vortex Filament from startpos to endpos with strength circulation. '''
    def __init__(self, startPos: Vec, endPos: Vec):
        self.startPos = startPos
        self.endPos = endPos
        self.centre = (endPos + startPos)/2

        self.circulation = None

    def set_circulation(self, magnitude):
        self.circulation = magnitude

    def GetInducedFlow(self, pos: Vec):
        '''Gets the induced flow due to this filament at the global position described by pos. Returns a Vec for the flow (self-centered).'''
        # Could be done with nice vectors, but just using the algorithm given by the slides
        d1 = pos - self.startPos
        d2 = pos - self.endPos
        d0 = self.endPos - self.startPos
        r1 = d1.Length()
        r2 = d2.Length()
        r12x = d1.yglob*d2.zglob - d1.zglob*d2.yglob
        r12y = -d1.xglob*d2.zglob + d1.zglob*d2.xglob
        r12z = d1.xglob*d2.yglob - d1.yglob*d2.xglob
        r12sqr = r12x*r12x + r12y+r12y + r12z*r12z
        r01 = d0.xglob*d1.xglob + d0.yglob*d1.yglob + d0.zglob*d1.zglob
        r02 = d0.xglob*d2.xglob + d0.yglob*d2.yglob + d0.zglob*d2.zglob

        K = self.circulation / (4 * np.pi * r12sqr) * (r01/r1 - r02/r2)

        return Vec((K*r12x, K*r12y, K*r12z))

    def plot(self, colour, ax=None):
        startPosXYZ = self.startPos.output_xyz(glob=True)
        endPosXYZ = self.endPos.output_xyz(glob=True)

        x, y, z = (startPosXYZ[0], endPosXYZ[0]), (startPosXYZ[1], endPosXYZ[1]), (startPosXYZ[2], endPosXYZ[2])
        if ax is None:
            plt.scatter(*startPosXYZ, color=colour)

        else:
            ax.scatter(*startPosXYZ, color=colour)
            ax.plot(x, y, z, color=colour)


class Leg:
    def __init__(self, reset=False):
        if not reset:
            self.control_points = []

        self.circulation = None

    def CreateControlPoints(self, arrayOfCylindricalPositions):
        '''Takes an argument like arrayOfCylindricalPositions = [(2,pi,0), (2.1,pi,1),...] and a root reference position in carthesian.'''
        
        # Requires a starting point!
        # TODO: Aida: Make this iterate on index and not in this stupid way with popping
        prevPoint = arrayOfCylindricalPositions[0]
        arrayOfCylindricalPositions.pop(0)
        for cylPos in arrayOfCylindricalPositions:
            self.control_points.append(Filament(cylPos, prevPoint))  
            prevPoint = cylPos     

    def reset(self):
        [cp.reset() for cp in self.control_points]
        self.__init__(reset=True)

    def plot(self, colour, ax=None):
        for filament in self.control_points:
            filament.plot(colour, ax=ax)


class HorseShoe:

    def __init__(self, airfoil, chord, pos_inner, pos_outer, twist, reset=False):
        if not reset:
            # TODO: Which leg is the outer one and which one is the inner one? Rename.
            self.leg_outer = Leg()
            self.leg_inner = Leg()

            self.airfoil = airfoil
            self.chord = chord
            self.twist = twist

            self.delta_r = (pos_outer - pos_inner).Length()
            self.pos_inner = pos_inner
            self.pos_outer = pos_outer
            self.pos_centre = (pos_inner + pos_outer) / 2  # vorticity needs to be calculated at the blade element center

        self.circulation = None
        self.induced_velocity = None

        self.vind_z = None
        self.vind_theta = None
        self.w_flow = None
        self.w_rot = None
        self.w = None
        self.phi = None
        self.alpha = None

    def set_circulation(self, v_inf, omega):
        self.vind_z = self.induced_velocity.zglob
        self.vind_theta = self.induced_velocity[1] * np.cos(self.pos_centre.thetaloc) - self.induced_velocity.xglob * np.sin(self.pos_centre.thetaloc)

        self.w_flow = v_inf + self.vind_z
        self.w_rot = omega * self.pos_centre.rloc + self.vind_theta

        self.w = np.sqrt(self.w_flow*self.w_flow + self.w_rot*self.w_rot)
        self.phi = np.arctan2(self.w_flow, self.w_rot)
        self.alpha = self.phi - self.twist

        self.circulation = .5 * self.w * self.delta_r * self.chord * self.airfoil.cl(np.degrees(self.alpha))

    def GetInducedVelocityInducedByHorseshoe(self, pos: Vec):
        '''Gets the total induced by the horseshoe at a specific point in space.'''
        totalInducedVelocity = Vec((0,0,0))

        for filament in self.leg_inner.control_points:
            flowByFilament = filament.GetInducedFlow(pos)
            totalInducedVelocity += flowByFilament

        for filament in self.leg_outer.control_points:
            flowByFilament = filament.GetInducedFlow(pos)
            totalInducedVelocity += flowByFilament

        return totalInducedVelocity

    def reset(self):
        self.leg_inner.reset()
        self.leg_outer.reset()
        self.__init__(0, 0, 0, 0, 0, reset=True)

    def plot(self, base_colour, ax=None):
        colour = base_colour * (self.pos_centre.rloc + 50) / (2*50)
        self.leg_inner.plot(colour, ax=ax)
        self.leg_outer.plot(colour, ax=ax)


    def GenerateWakePoints(self, axialVelocity: float, rotSpeed: float, resolution: int, tmax: float):
        '''Takes in an axial velocity after the turbine to consider previous points the turbine would have been in in the past.'''

        outer_leg_control_points = [];
        inner_leg_control_points = [];

        intervals = np.linspace(0, tmax, resolution)
        for t in intervals:
            # Skip the 0 time, because that's where stuff's currently at?
            # TODO: check if the actual vortices generated by the airfoil at dt = 0 are accounted for.

            z = t*axialVelocity

            # Every horseshoe has a pos_inner and pos_outer, so do the calculation twice
            y_inner = self.pos_inner.rloc * np.sin(rotSpeed * t + self.pos_inner.thetaloc)
            x_inner = self.pos_inner.rloc * np.cos(rotSpeed * t + self.pos_inner.thetaloc)
            pos_inner_trail = Vec((x_inner, y_inner, z))
            inner_leg_control_points.append(pos_inner_trail)

            y_outer = self.pos_outer.rloc * np.sin(rotSpeed * t + self.pos_outer.thetaloc)
            x_outer = self.pos_outer.rloc * np.cos(rotSpeed * t + self.pos_outer.thetaloc)
            pos_outer_trail = Vec((x_outer, y_outer, z))
            outer_leg_control_points.append(pos_outer_trail)
        
        # With the points generated, hand over to each leg to convert the control points into filaments
        self.leg_outer.CreateControlPoints(outer_leg_control_points)
        self.leg_inner.CreateControlPoints(inner_leg_control_points)


            

class Turbine:
    """
    Turbine parameters, air density, U_inf
    """
    def __init__(self, rotation=0, referencePos=(0.,0.,0.), reset=False):
        # data = read_from_file('DU95W150.csv')
        # self.alpha_lst = data[:, 0]
        # self.cl_lst = data[:, 1] #; self.cd_lst = data[:, 2]; self.cm_lst = data[:, 3]

        self.b = 3 # Number of blades
        self.wakePointResolution = 10
        self.twakemax = 2.5
        self.n_elements = 3 # Divide the blade up in n_elements
        self.rho = 1.225 # [Kg/m3]
        self.u_inf = 10 # [m/s] U_infinity = free stream velocity
        self.radius = 50 # Total radius
        self.tsr = 10
        self.omega = self.tsr * self.u_inf / self.radius
        self.blade_pitch = 0
        r_start = 0.2*self.radius

        airfoil = DU95W150()

        self.rotation = rotation

        self.horseshoes = [[], [], []]

        for i in range(self.n_elements):
            r_inner = r_start + (self.radius - r_start) / self.n_elements * i
            r_outer = r_start + (self.radius - r_start) / self.n_elements * (i+1)
            r = np.mean([r_inner, r_outer])
            # Sorry for hardcoding the equations below- taken from the assignment description :)
            twist = 14 * (1 - r / self.radius)
            chord = (3 * (1 - r / self.radius) + 1)

            # BladeElement takes in argument relative_pitch, I assume that this means total? So offset with the blade pitch
            relative_pitch = self.blade_pitch + twist

            # Horseshoes is a 2d array of horseshoes, one for each rotor. This iterates over each rotor.
            for idx, _ in enumerate(self.horseshoes):
                pos_inner = Vec((r_inner, rotation + idx * np.radians(120), 0), referencePos=referencePos, bLocalCylindrical=True)
                pos_outer = Vec((r_outer, rotation + idx * np.radians(120), 0), referencePos=referencePos, bLocalCylindrical=True)

                horseshoeToAdd = HorseShoe(airfoil, chord, pos_inner, pos_outer, relative_pitch)
                horseshoeToAdd.GenerateWakePoints(self.u_inf, self.omega, self.wakePointResolution, self.twakemax)
                
                self.horseshoes[idx].append(horseshoeToAdd)

    def GetInducedVelocityByTurbine(self, pos: Vec):
        '''Iterates over every horseshoe, gets their induced velocity and plunges them all together.'''
        totalInducedVelocity = Vec((0,0,0))

        for horseshoe in self.horseshoes:
            inducedVelocityByHorseshoe = horseshoe.GetInducedVelocityInducedByHorseshoe(pos)
            totalInducedVelocity += inducedVelocityByHorseshoe

        return totalInducedVelocity

    def reset(self):
        [[hs.reset() for hs in blade] for blade in self.horseshoes]
        self.__init__(self.rotation, reset=True)

    def SetInducedVelocityForHorseshoes(self):
        '''For each horseshoe, sets the induced velocity by all the other horseshoes and itself at its centrepoint.'''
        # Iterate over the horseshoes to set
        for set in self.horseshoes:
            set.induced_velocity = self.GetInducedVelocityByTurbine(set.pos_centre)

    def set_circulations_horseshoes(self):
        for blade in self.horseshoes:
            for set in blade:
                set.set_circulation(self.u_inf, self.omega)

    def plot(self, ax=None):
        for i, horseshoes in enumerate(self.horseshoes):
            colour = np.zeros(3)
            colour[i] = 1
            for horseshoe in horseshoes:
                horseshoe.plot(colour, ax=ax)


if __name__ == "__main__":
    print("This is a lifting line library, pls dont run this")
    #Turbine()
