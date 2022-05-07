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
    def __init__(self, startPos, endPos):
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


class Leg:
    def __init__(self, reset=False):
        if not reset:
            self.control_points = []

        self.circulation = None

    def CreateControlPoints(self, arrayOfCylindricalPositions: list(Vec), rootPos):
        '''Takes an argument like arrayOfCylindricalPositions = [(2,pi,0), (2.1,pi,1),...] and a root reference position in carthesian.'''
        
        # Requires a starting point!
        # TODO: Aida: Make this iterate on index and not in this stupid way with popping
        prevPoint = arrayOfCylindricalPositions[0]
        arrayOfCylindricalPositions.pop(0)
        for cylPos in arrayOfCylindricalPositions:
            point = Vec((cylPos.x, cylPos.y, cylPos.z), rootPos, bLocalCylindrical=True)
            self.control_points.append(Filament(point, prevPoint))  
            prevPoint = point     

    def reset(self):
        [cp.reset() for cp in self.control_points]
        self.__init__(reset=True)


class HorseShoe:

    def __init__(self, airfoil, chord, pos_inner, pos_outer, twist, reset=False):
        if not reset:
            # TODO: Which leg is the outer one and which one is the inner one? Rename.
            self.leg1 = Leg()
            self.leg2 = Leg()

            self.airfoil = airfoil
            self.chord = chord
            self.twist = twist

            self.delta_r = (pos_outer - pos_inner).Length()
            self.pos_inner = pos_inner
            self.pos_outer = pos_outer
            self.pos_centre = (pos_inner + pos_outer) / 2  # vorticity needs to be calculated at the blade element center

        self.circulation = None
        self.induced_velocity = None

    def set_circulation(self, v_inf, omega):
        vind_z = self.induced_velocity.zglob
        vind_theta = self.induced_velocity[1] * np.cos(self.pos_centre.thetaloc) - self.induced_velocity.xglob * np.sin(self.pos_centre.thetaloc)

        w_flow = v_inf + vind_z
        w_rot = omega * self.pos_centre.rloc + vind_theta

        w = np.sqrt(w_flow*w_flow + w_rot*w_rot)
        phi = np.arctan2(w_flow, w_rot)

        self.circulation = .5 * w * self.delta_r * self.chord * self.airfoil.cl(np.degrees(phi - self.twist))

    def GetInducedVelocityInducedByHorseshoe(self, pos: Vec):
        '''Gets the total induced by the horseshoe at a specific point in space.'''
        totalInducedVelocity = Vec((0,0,0))

        for filament in self.leg1.control_points:
            flowByFilament = filament.GetInducedFlow(pos)
            totalInducedVelocity += flowByFilament

        for filament in self.leg2.control_points:
            flowByFilament = filament.GetInducedFlow(pos)
            totalInducedVelocity += flowByFilament

        return totalInducedVelocity

    def reset(self):
        self.leg1.reset()
        self.leg2.reset()
        self.__init__(0, 0, 0, 0, 0, reset=True)

    def GenerateWakePoints(self, axialVelocity: float, rotSpeed: float, resolution: int, tmax: float):
        '''Takes in an axial velocity after the turbine to consider previous points the turbine would have been in in the past.'''

        outer_leg_control_points = [];
        inner_leg_control_points = [];

        intervals = np.linspace(0, tmax, resolution)
        for t in intervals:
            # Skip the 0 time, because that's where stuff's currently at?
            # TODO: check if the actual vortices generated by the airfoil at dt = 0 are accounted for.
            if t == 0:
                continue;

            dx = t*axialVelocity

            # Every horseshoe has a pos_inner and pos_outer, so do the calculation twice
            dy_inner = self.pos_inner.r * np.sin(rotSpeed * t)
            dz_inner = self.pos_inner.r * np.cos(rotSpeed * t)
            pos_inner_trail = self.pos_inner + Vec((dx, dy_inner, dz_inner))
            inner_leg_control_points.append(pos_inner_trail)

            dy_outer = self.pos_outer.r * np.sin(rotSpeed * t)
            dz_outer = self.pos_outer.r * np.cos(rotSpeed * t)
            pos_outer_trail = self.pos_outer + Vec((dx,dy_outer, dz_outer))
            outer_leg_control_points.append(pos_outer_trail)
        
        # With the points generated, hand over to each leg to convert the control points into filaments
        #TODO: Check if these feed into the right legs.
        self.leg1.CreateControlPoints(outer_leg_control_points)
        self.leg2.CreateControlPoints(inner_leg_control_points)


            

class Turbine:
    """
    Turbine parameters, air density, U_inf
    """
    def __init__(self, rotation, reset=False):
        # data = read_from_file('DU95W150.csv')
        # self.alpha_lst = data[:, 0]
        # self.cl_lst = data[:, 1] #; self.cd_lst = data[:, 2]; self.cm_lst = data[:, 3]

        self.b = 3 # Number of blades
        self.n_elements = 1 # Divide the blade up in n_elements
        self.rho = 1.225 # [Kg/m3]
        self.u_inf = 10 # [m/s] U_infinity = free stream velocity
        self.radius = 50 # Total radius
        self.tsr = 10
        self.omega = self.tsr * self.u_inf / self.radius
        self.blade_pitch = 0
        r_start = 0.2*self.radius

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

            for idx, _ in enumerate(self.horseshoes):
                pos_inner = Vec((r_inner, rotation + idx * np.radians(120), 0))
                pos_outer = Vec((r_outer, rotation + idx * np.radians(120), 0))
                self.horseshoes[idx].append(HorseShoe(DU95W150, chord, pos_inner, pos_outer, relative_pitch))

    def GetInducedVelocityByTurbine(self, pos: Vec):
        '''Iterates over every horseshoe, gets their induced velocity and plunges them all together.'''
        totalInducedVelocity = Vec((0,0,0))

        for horseshoe in self.horseshoes:
            inducedVelocityByHorseshoe = horseshoe.GetInducedVelocityInducedByHorseshoe(pos)
            totalInducedVelocity += inducedVelocityByHorseshoe

        return totalInducedVelocity

    def reset(self):
        [hs.reset() for hs in self.horseshoes]
        self.__init__(reset=True)

    def SetInducedVelocityForHorseshoes(self):
        '''For each horseshoe, sets the induced velocity by all the other horseshoes and itself at its centrepoint.'''
        # Iterate over the horseshoes to set
        for set in self.horseshoes:
            set.induced_velocity = self.GetInducedVelocityByTurbine(set.pos_centre)

    def set_circulations_horseshoes(self):
        for blade in self.horseshoes:
            for set in blade:
                set.set_circulation(self.u_inf, self.omega)


if __name__ == "__main__":
    print("This is a lifting line library, pls dont run this")
    Turbine()
