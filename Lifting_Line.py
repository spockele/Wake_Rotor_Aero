from logging import root
from multiprocessing.sharedctypes import Value
import numpy as np
from BEM_code import DU95W150
import matplotlib.pyplot as plt
from read_write import read_from_file, write_to_file
import scipy.integrate as spig


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
            raise TypeError("b-baka, are you sure you passed a tuple into the constructor for Vec like Vec((1,2,3)), and not Vec(1,2,3)?")

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
        # Ugly copy all the thingssss
        xp = pos.xglob
        yp = pos.yglob
        zp = pos.zglob

        x1 = self.startPos.xglob
        y1 = self.startPos.yglob
        z1 = self.startPos.zglob

        x2 = self.endPos.xglob
        y2 = self.endPos.yglob
        z2 = self.endPos.zglob

        r1 = np.sqrt((xp - x1)**2 + (yp - y1)**2 + (zp - z1)**2)
        r2 = np.sqrt((xp - x2)**2 + (yp - y2)**2 + (zp - z2)**2)
        r12x = (yp - y1)*(zp - z2) - (zp - z1)*(yp - y2)
        r12y = -1*(xp - x1)*(zp - z2) + (zp - z1)*(xp - x2)
        r12z = (xp - x1)*(yp - y2) - (yp - y1)*(xp - x2)
        r12sqr = r12x**2 + r12y**2 + r12z**2
        r01 = (x2 - x1)*(xp - x1) + (y2 - y1)*(yp - y1) + (z2 - z1)*(zp - z1)
        r02 = (x2 - x1)*(xp - x2) + (y2 - y1)*(yp - y2) + (z2 - z1)*(zp - z2)
        if abs(r12sqr) < 1e-9:
            return Vec((0,0,0))
        k = self.circulation/(4*np.pi*r12sqr)*(r01/r1 - r02/r2)
        inducedVel = Vec((k*r12x, k*r12y, k*r12z))
        return inducedVel


    def plot(self, colour, ax=None):
        startPosXYZ = self.startPos.output_xyz(glob=True)
        endPosXYZ = self.endPos.output_xyz(glob=True)

        x, y, z = (startPosXYZ[0], endPosXYZ[0]), (startPosXYZ[1], endPosXYZ[1]), (startPosXYZ[2], endPosXYZ[2])
        if ax is None:
            plt.scatter(*startPosXYZ, marker='.', color=colour)

        else:
            # ax.scatter(*startPosXYZ, marker='.', color=colour)
            ax.plot(z, x, y, color=colour)


class Leg:
    def __init__(self):
        self.control_points = []
        self.circulation = None

    def CreateControlPoints(self, arrayOfCylindricalPositions):
        '''Takes an argument like arrayOfCylindricalPositions = [(2,pi,0), (2.1,pi,1),...] and a root reference position in carthesian.'''
        
        # Requires a starting point!
        # TODO: Aida: Make this iterate on index and not in this stupid way with popping
        prevPoint = arrayOfCylindricalPositions[0]
        for cylPos in arrayOfCylindricalPositions[1:]:
            self.control_points.append(Filament(prevPoint, cylPos))
            prevPoint = cylPos

    def plot(self, colour, ax=None):
        for filament in self.control_points:
            filament.plot(colour, ax=ax)


class HorseShoe:

    def __init__(self, airfoil, chord, pos_inner, pos_outer, twist):
        self.leg_outer = Leg()
        self.leg_inner = Leg()

        self.airfoil = airfoil
        self.chord = chord
        self.twist = twist

        self.delta_r = (pos_outer - pos_inner).Length()
        self.pos_inner = pos_inner
        self.pos_outer = pos_outer
        self.pos_centre = (pos_inner + pos_outer) / 2  # vorticity needs to be calculated at the blade element center

        # Filament that actually represents the lifting line
        self.liftingLine = Filament(self.pos_inner, self.pos_outer)

        self.circulation = None
        self.induced_velocity = Vec((0,0,0))

        self.vind_z = None
        self.vind_theta = None
        self.w_flow = None
        self.w_rot = None
        self.w = None
        self.phi = None
        self.alpha = None
        self.pt = None
        self.pn = None
        self.a = None # Axial induction factor
        self.a_prime = None # Tangential induction factor

    def set_circulation(self, v_inf, omega, relaxationFactor):
        self.vind_z = self.induced_velocity.zglob
        self.vind_theta = self.induced_velocity.yglob * np.cos(self.pos_centre.thetaloc) - self.induced_velocity.xglob * np.sin(self.pos_centre.thetaloc)

        self.w_flow = v_inf + self.vind_z

        # @fien: see what making this plus a minus does
        self.w_rot = omega * self.pos_centre.rloc + self.vind_theta

        self.w = np.sqrt(self.w_flow*self.w_flow + self.w_rot*self.w_rot)
        self.phi = np.arctan2(self.w_flow, self.w_rot)
        self.alpha = self.phi - self.twist

        # Little hack to make sure that the first pre-iteration doesn't return all nulls
        if self.circulation == None:
            previousCirculation = 0
        else:
            previousCirculation = self.circulation

        # Since static iteration tends to wildly oscillate, we add a relaxation factor.
        # We make this relaxation factor a function of the gradient of the circulation (first order newton)
        self.circulation = relaxationFactor * (.5 * self.w * self.delta_r * self.chord * self.airfoil.cl(np.degrees(self.alpha))) + (1-relaxationFactor) * previousCirculation

        # Propogate the circulation over to all the filaments in this horseshoe
        # TODO: ABSOLUTELY DOUBLE CHECK IF WE MUTLIPLY THE RIGHT ONE WITH -1
        self.liftingLine.set_circulation(self.circulation)

        for filament in self.leg_inner.control_points:
            filament.set_circulation(self.circulation)
        for filament in self.leg_outer.control_points:
            filament.set_circulation(-1*self.circulation)

        return self.circulation - previousCirculation, self.delta_r

    def GetForcesAndFactors(self, rho, v_inf, omega):
        self.LiftDensity = 0.5*rho*self.w*self.w* self.chord*self.airfoil.cl(np.degrees(self.alpha))
        self.DragDensity = 0.5*rho*self.w*self.w* self.chord*self.airfoil.cd(np.degrees(self.alpha))
        self.pn = self.DragDensity*np.sin(self.phi)+self.LiftDensity*np.cos(self.phi)
        self.pt = self.LiftDensity*np.sin(self.phi)-self.DragDensity*np.cos(self.phi)
        self.a = 1-self.w_flow/v_inf
        self.a_prime = self.w_rot/omega/self.pos_centre.rloc - 1

    def GetInducedVelocityInducedByHorseshoe(self, pos: Vec):
        '''Gets the total induced by the horseshoe at a specific point in space.'''
        totalInducedVelocity = Vec((0,0,0))

        # Do NOT forget the contribution by the actual lifting line lmao
        totalInducedVelocity += self.liftingLine.GetInducedFlow(pos)

        for filament in self.leg_inner.control_points:
            flowByFilament = filament.GetInducedFlow(pos)
            totalInducedVelocity += flowByFilament

        for filament in self.leg_outer.control_points:
            flowByFilament = filament.GetInducedFlow(pos)
            totalInducedVelocity += flowByFilament

        return totalInducedVelocity

    def plot(self, base_colour, ax=None):
        colour = base_colour * (self.pos_centre.rloc + 50) / (2*50)
        self.leg_inner.plot(colour, ax=ax)
        self.leg_outer.plot(colour, ax=ax)


    def GenerateWakePoints(self, axialVelocity: float, rotSpeed: float, resolution: int, tmax: float):
        '''Takes in an axial velocity after the turbine to consider previous points the turbine would have been in in the past.'''

        outer_leg_control_points = []
        inner_leg_control_points = []

        intervals = np.linspace(0, tmax, resolution)
        for t in intervals:
            # Skip the 0 time, because that's where stuff's currently at?
            # TODO: check if the actual vortices generated by the airfoil at dt = 0 are accounted for.

            if not t:
                inner_leg_control_points.append(self.pos_inner)
                outer_leg_control_points.append(self.pos_outer)

            z = t*axialVelocity + self.pos_inner.refPos[2] + self.chord * np.sin(self.twist)
            az_offset = self.chord * np.cos(self.twist)
            r_inner = np.sqrt(self.pos_inner.rloc*self.pos_inner.rloc + az_offset*az_offset)
            r_outer = np.sqrt(self.pos_outer.rloc*self.pos_outer.rloc + az_offset*az_offset)

            theta_offset_inner = np.arctan2(az_offset, self.pos_inner.rloc)
            theta_offset_outer = np.arctan2(az_offset, self.pos_outer.rloc)

            # Every horseshoe has a pos_inner and pos_outer, so do the calculation twice
            y_inner = r_inner * np.sin(rotSpeed * t + self.pos_inner.thetaloc + theta_offset_inner) + self.pos_inner.refPos[1]
            x_inner = r_inner * np.cos(rotSpeed * t + self.pos_inner.thetaloc + theta_offset_inner) + self.pos_inner.refPos[0]
            pos_inner_trail = Vec((x_inner, y_inner, z))
            inner_leg_control_points.append(pos_inner_trail)

            y_outer = r_outer * np.sin(rotSpeed * t + self.pos_outer.thetaloc + theta_offset_outer) + self.pos_inner.refPos[1]
            x_outer = r_outer * np.cos(rotSpeed * t + self.pos_outer.thetaloc + theta_offset_outer) + self.pos_inner.refPos[0]
            pos_outer_trail = Vec((x_outer, y_outer, z))
            outer_leg_control_points.append(pos_outer_trail)
        
        # With the points generated, hand over to each leg to convert the control points into filaments
        self.leg_outer.CreateControlPoints(outer_leg_control_points)
        self.leg_inner.CreateControlPoints(inner_leg_control_points)


class Turbine:
    """
    Turbine parameters, air density, U_inf
    """
    def __init__(self, n_rot_wake=8, n_point_per_rotation=12, n_blade_elements=30, convection_speed=10, tsr=10, rotation=0, referencePos=(0.,0.,0.)):
        self.b = 3 # Number of blades
        self.radius = 50  # Total radius
        self.blade_pitch = np.radians(-2)
        self.rotation = rotation

        self.rho = 1.225 # [Kg/m3]
        self.u_inf = 10 # [m/s] U_infinity = free stream velocity
        self.tsr = tsr

        self.omega = self.tsr * self.u_inf / self.radius
        r_start = 0.2*self.radius
        airfoil = DU95W150()

        self.wakePointResolution = n_rot_wake * n_point_per_rotation
        self.twakemax = n_rot_wake * 2 * np.pi / self.omega
        self.n_elements = n_blade_elements  # Divide the blade up in n_elements
        self.u_wake = convection_speed

        self.horseshoes = [[], [], []]
        for i in range(self.n_elements):
            r_inner = r_start + (self.radius - r_start) / self.n_elements * i
            r_outer = r_start + (self.radius - r_start) / self.n_elements * (i+1)

            r = np.mean([r_inner, r_outer])
            # Sorry for hardcoding the equations below- taken from the assignment description :)
            twist = np.radians(14 * (1 - r / self.radius))
            chord = (3 * (1 - r / self.radius) + 1)

            # BladeElement takes in argument relative_pitch, I assume that this means total? So offset with the blade pitch
            relative_pitch = self.blade_pitch + twist

            # Horseshoes is a 2d array of horseshoes, one for each rotor. This iterates over each rotor.
            for idx, _ in enumerate(self.horseshoes):
                pos_inner = Vec((r_inner, rotation + idx * np.radians(120), 0), referencePos=referencePos, bLocalCylindrical=True)
                pos_outer = Vec((r_outer, rotation + idx * np.radians(120), 0), referencePos=referencePos, bLocalCylindrical=True)

                horseshoeToAdd = HorseShoe(airfoil, chord, pos_inner, pos_outer, relative_pitch)
                horseshoeToAdd.GenerateWakePoints(self.u_wake, self.omega, self.wakePointResolution, self.twakemax)
                
                self.horseshoes[idx].append(horseshoeToAdd)

    def GetInducedVelocityByTurbine(self, pos: Vec):
        '''Iterates over every horseshoe, gets their induced velocity and plunges them all together.'''
        totalInducedVelocity = Vec((0,0,0))

        for rotor in self.horseshoes:
            for horseshoe in rotor:
                inducedVelocityByHorseshoe = horseshoe.GetInducedVelocityInducedByHorseshoe(pos)
                totalInducedVelocity += inducedVelocityByHorseshoe

        return totalInducedVelocity

    def SetInducedVelocityForHorseshoes(self, turbines=None):
        '''For each horseshoe, sets the induced velocity by all the other horseshoes and itself at its centrepoint.'''
        # Iterate over the horseshoes to set
        if turbines is None:
            turbines = [self]

        for blade in self.horseshoes:
            for horseshoe in blade:
                new_induced_velocity = Vec((0, 0, 0))
                for turbine in turbines:
                    new_induced_velocity += turbine.GetInducedVelocityByTurbine(horseshoe.pos_centre)

                horseshoe.induced_velocity = new_induced_velocity

    def set_circulations_horseshoes(self, relaxationFactor):
        '''Sets the circulation for all the horseshoes based on their internally saved flow deviation vector. Returns the change in circulation (delta gamma) for the element that has it as the highest.'''
        highestDeltaGamma = 0
        highestIndex = 0, 0
        dr = 0
        for j, blade in enumerate(self.horseshoes):
            for i, set in enumerate(blade):
                # Calls function that updates, and returns the change.
                deltaGamma, dr = set.set_circulation(self.u_inf, self.omega, relaxationFactor)
                if abs(deltaGamma) > abs(highestDeltaGamma):
                    highestDeltaGamma = deltaGamma
                    dr = dr
                    highestIndex = j, i

        return highestIndex, highestDeltaGamma, dr

    def plot(self, ax=None):
        for i, horseshoes in enumerate(self.horseshoes):
            colour = np.zeros(3)
            colour[i] = 1
            for horseshoe in horseshoes:
                horseshoe.plot(colour, ax=ax)

    def extract_information_N_write(self, write=True, suffix=''):
        out_array = np.empty((3, 7, len(self.horseshoes[0])))
        thrust = 0
        power = 0
        for i, blade in enumerate(self.horseshoes):
            for j, horseshoe in enumerate(blade):
                horseshoe.GetForcesAndFactors(self.rho, self.u_inf, self.omega)
                out_array[i, 0, j] = horseshoe.pos_centre.rloc
                out_array[i, 1, j] = horseshoe.alpha
                out_array[i, 2, j] = horseshoe.phi
                out_array[i, 3, j] = horseshoe.pn
                out_array[i, 4, j] = horseshoe.pt
                out_array[i, 5, j] = horseshoe.a
                out_array[i, 6, j] = horseshoe.a_prime
            if write:
                write_to_file(out_array[i,:,:], f'saved_data/LL_output_{self.tsr}_blade{i}{suffix}.txt')
            # Thrust and Power calculation
            thrust += spig.trapz(out_array[i, 3,:], out_array[i, 0,:])
            power += self.omega * spig.trapz(out_array[i, 4,:] * out_array[i, 0,:], out_array[i, 0,:])

        CT = thrust / (0.5 * self.rho * np.pi * self.radius ** 2 * self.u_inf ** 2)
        Cp = power / (0.5 * self.rho * np.pi * self.radius ** 2 * self.u_inf ** 3)

        if write:
            write_to_file([[CT, Cp]], f'saved_data/LL_cT_cp_{self.tsr}{suffix}.txt')

        return out_array, CT, Cp


def GetRelaxationFactor(deltaCirculation, dr):
    # gamma_n = f * gamma + (1-f)*gamma_{n-1}

    # 1 |
    #   |-K
    #   |  \
    #   |   L---
    # 0 |_______

    ky = 0.2 # K_y
    kx = .11 / dr # K_x
    ly = 0.5 # L_y
    lx = 25 / dr # L_x

    if deltaCirculation is None:
        return ky
    elif deltaCirculation < kx:
        f = ky
    elif deltaCirculation > lx:
        f = ly
    else:
        slope = (ky - ly)/(kx - lx)
        offset = ky - slope * kx
        f = slope * abs(deltaCirculation) + offset
    return f


if __name__ == "__main__":
    print("This is a lifting line library, pls dont run this")
