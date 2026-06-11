import math
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as LA


#Test forces with simple chain
#Compute manually as well to compare

#Dynamics


class Bead():
    def __init__(self, mass, hydro, *args):
        '''
        Bead initializer
        mass -> bead mass (relative)
        hydro -> char, Hydro Neutral(N), Phillic (L), Phobi (B)
        args:
        Either coordinates as 3 doubles
        Or 1 array of 3 coordinates
        Or 2 arrays, 1 for coordinate, 1 for velocity
        '''

        self.mass=mass
        self.hydro=hydro

        self.velocity = np.zeros(3)

        if len(args)==3:
            self.coords = np.array([args[0], args[1], args[2]])

        if len(args)==1:
            self.coords = np.array(args[0])

        if len(args)==2:
            self.coords = np.array(args[0])
            self.velocity = np.array(args[1])

    def __str__(self):
        return f" mass = {self.mass} \n hydro {self.hydro} \n coords ({self.coords[0]},{self.coords[1]},{self.coords[2]})"

    def __getitem__(self, idx):
        return self.coords[idx]

    def __setitem__(self, idx, value):
        self.coords[idx]=value

    def get_momentum(self):
        return self.mass*self.velocity

    def plot(self, ax, resolution=32, bead_scale=0.1, **kwargs):
        radius = self.mass

        color_map = {
            "N": "green",
            "B": "red",
            "L": "blue"
        }

        u = np.linspace(0,2*np.pi, resolution)
        v = np.linspace(0, np.pi, resolution)

        x = self.coords[0] + (radius + np.outer(np.cos(u), np.sin(v)))*bead_scale
        y = self.coords[1] + (radius + np.outer(np.sin(u), np.sin(v)))*bead_scale
        z = self.coords[2] + (radius + np.outer(np.ones_like(u), np.cos(v)))*bead_scale

        ax.plot_surface(x,y,z, color=color_map.get(self.hydro, "gray"),
                        **kwargs)

def angle_between(v1, v2):
    #Returns the acute angles between to 3D vectors relative to the normal of their plane
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    return abs(np.arctan2(np.linalg.norm(np.cross(v1_u, v2_u)), np.dot(v1_u, v2_u)))

class Chain():
    def __init__(self, *args):
        '''
        A Chain is a set of beads connected in order (as related by indexes)
        The initializer either generates an empty array, or generate a chain from
        an array of beads

        The length of the chain is N=len(self.beads)

        chain_vector index k is the chain vector  beads[k+1]-beads[k], size N-1

        The bond angle k is the acute angle between the coordinates of beads
        k,k+1,k+2

        The dihedral_angle k is the dihedral torsion angle defined by the coordinates of
        beads k,k+1,k+2,k+3

        if a method updates the chain or adds a bead, the flag update_chain is set to true


        '''
        self.chain_vector = []

        self.bond_angles = []
        self.dihedral_angles = []
        self.update_chain=True
        if len(args)>0:
            self.beads = args[0]
            self.update_chain_vec()
        else:
            self.beads = []

    def update_chain_vec(self, force=False):
        '''
        This update the chain with chain vectors, bond angle and dihedral angle
        '''
        if self.update_chain==True or force:
            self.update_chain=False
            BL_idx=0
            BA_idx=0
            DA_idx=0
            for idx in range(1, len(self.beads)):
                r1 = self.beads[idx-1].coords
                r2 = self.beads[idx].coords

                #Chain vectors are used over and over in differen computations, so we keep track of them
                u_current = r2-r1

                if idx < len(self.chain_vector):
                    self.chain_vector[idx] = u_current
                else:
                    self.chain_vector.append(u_current)

                BL_idx+=1

                #Get Bond Angles
                if idx>1:
                    bond_angle = angle_between(-u_prev, u_current)
                    if BA_idx < len(self.bond_angles):
                        self.bond_angles[BA_idx] = bond_angle
                    else:
                        self.bond_angles.append(bond_angle)
                    BA_idx+=1

                #Get Dihedral angles
                if idx>2:
                    n_0 = np.cross(-u_prev_prev, u_prev)
                    n_1 = np.cross(-u_prev, u_current)
                    DA_angle=0
                    if LA.norm(n_0)*LA.norm(n_1)>1e-16:
                        DA_angle=angle_between(n_0, n_1)

                    if DA_idx<len(self.dihedral_angles):
                        self.dihedral_angles[DA_idx] = DA_angle
                    else:
                        self.dihedral_angles.append(DA_angle)
                    DA_idx+=1

                if idx==1:
                    u_prev = u_current
                    u_prev_prev = u_current
                else:
                    u_prev_prev = u_prev
                    u_prev = u_current


    def __len__(self):
        return len(self.beads)

    def __getitem__(self, idx):
        return self.beads[idx]

    def __str__(self):
        chain_str=""
        for idx in range(len(self.beads)):
            chain_str+= f"Bead {idx} \n{self.beads[idx]} \n"

        return chain_str

    def plot(self, ax, resolution=32, bead_scale=0.1, **kwargs):
        for bead in self.beads:
            bead.plot(ax, resolution, bead_scale, **kwargs)

        for b1, b2 in zip(self.beads[:-1], self.beads[1:]):
            x = [b1[0], b2[0]]
            y = [b1[1], b2[1]]
            z = [b1[2], b2[2]]

            ax.plot(x,y,z, color="black", linewidth=2)
# ============================================================================= #
#                                POTENTIAL FUNCTIONS
# ============================================================================= #


class Hydro_Chain(Chain):
    def __init__(self, beads, hydro_strengh=1, avg_bond_length=1.2, lambda_std=1.0):
        super().__init__(beads)
        self.e_h=hydro_strengh
        self.a=avg_bond_length
        self.lambda_std=1.0
        self.forces = np.zeros((len(self.beads),3))

        self.bond_forces = self.forces
        self.bond_k = 100*self.e_h/(self.a*self.a)

        self.bond_angle_forces = self.forces
        self.angle_k = 20*self.e_h #per rad^2
        self.theta_0 = 1.8326 #rad

        self.dihedral_forces = self.forces


    def fill_bond_forces(self):

        for k in range(len(chain.chain_vector)):
            u = chain.chain_vector[k]
            u_abs = LA.norm(u)
            u_norm = u/u_abs

            f_mag = -2*self.bond_k*(u_abs - self.a)
            f_vec = u_norm*f_mag

            self.bond_forces[k] += f_vec
            self.bond_forces[k+1] += -f_vec

    def aux_get_corner_force(self, theta, f_dir):

        dir_norm = LA.norm(f_dir)

        if dir_norm < 1e-15:
            return [0,0,0]

        f_dir_norm = f_dir/dir_norm

        f = -2*self.angle_k*(theta - self.theta_0)

        return f_dir_norm*f


    def fill_angle_forces(self):

        for k in range(len(chain.bond_angles)):

            theta = self.bond_angles[k]
            u_0 = self.chain_vector[k]
            u_1 = self.chain_vector[k+1]

            u_0_x_u_1 = LA.cross(u_0, u_1)

            f_0_dir = LA.cross(u_0, u_0_x_u_1)
            f_1_dir = LA.cross(-u_1, u_0_x_u_1)

            f_0 = self.aux_get_corner_force(theta, f_0_dir)/LA.norm(u_0)
            f_1 = self.aux_get_corner_force(theta, f_1_dir)/LA.norm(u_1)

            self.bond_angle_forces[k] = f_0
            self.bond_angle_forces[k+1] = -f_1-f_0
            self.bond_angle_forces[k+2] = f_1


    def aux_dU_dtheta(self, A,B, theta):
        return 0.5*(A*np.sin(theta) + B*np.sin(3*theta))

    def fill_dihedral_forces(self):
        for k in range(len(chain.dihedral_angles)):

            BA_0 = self.bond_angles[k]
            BA_1 = self.bond_angles[k+1]

            sin_0 = np.sin(BA_0)
            sin_1 = np.sin(BA_1)

            if sin_0*sin_1 < 1e-15:
                return

            hydro_types = ""
            for jdx in range(4):
                hydro_types+=chain[k+jdx].hydro

            A = 1.2*self.e_h
            B=A

            if hydro_types.count("N") >=2:
                A=0
                B=0.2*self.e_h


            theta = self.dihedral_angles[k]

            u_0 = self.chain_vector[k]
            u_1 = self.chain_vector[k+1]
            u_2 = self.chain_vector[k+2]


            u_1Xu_0 = LA.cross(u_1, u_0)
            f_0_dir = u_1Xu_0/LA.norm(u_1Xu_0)
            f_0 = f_0_dir*self.aux_dU_dtheta(A,B, theta)/(LA.norm(u_0))

            u_1Xu_2 = LA.cross(u_1, u_2)
            f_3_dir = u_1Xu_2/LA.norm(u_1Xu_2)
            f_3 = f_3_dir*self.aux_dU_dtheta(A,B, theta)/(LA.norm(u_2))

            r_2 = 0.5*u_1
            tau_2 = -(LA.cross(r_2 +0.5*u_2, f_3) - 0.5*LA.cross(u_0, f_0))

            f_2 = LA.cross(tau_2, r_2)/(LA.norm(r_2)**2)

            self.dihedral_forces[k] = f_0
            self.dihedral_forces[k+1] = -f_0-f_2-f_3
            self.dihedral_forces[k+2] = f_2
            self.dihedral_forces[k+3] = f_3





def compute_V_DIH(chain, e_h=1):
    V_DIH = 0
    for idx in range(len(chain.chain_vector)-2):

        hydro_types = ""
        for jdx in range(4):ro_types+=chain[idx+jdx].hydro

        A=1.2*e_h
        B=A

        if hydro_types.count("N") >=2:
            A=0
            B=0.2*e_h

        dih_angle = chain.dihedral_angles[idx]
        cos_dihedral = np.cos(dih_angle)
        cos_3dihedral=np.cos(3*dih_angle)

        V_DIH += A*(1+cos_dihedral) + B*(1+cos_3dihedral)


    return V_DIH


def compute_V_BA(chain, e_h=1):
    k = 10*e_h
    theta_0 = 1.8326

    V_BA = 0

    for theta in chain.bond_angles:
        V_BA+= np.square(theta - theta_0)

    return V_BA

def compute_V_BL(chain, a=1 , e_h=1):
    k = 50*e_h/(a*a)
    V_BL=0
    for u in chain.chain_vector:
        V_BL+= k*np.square((LA.norm(u) - a))

    return V_BL


def compute_V_non(chain, a=1, e_h=1, lambda_std=1):

    e_L = (0.6667*e_h)
    N_beads = len(chain)
    V_non=0
    for idx in range(0, N_beads-3):
        for jdx in range(idx+3, N_beads):

            r = LA.norm(chain.beads[jdx].coords - chain.beads[idx].coords)

            if r > 1e-16:

                r_power_6 = np.power(a/r, 6)
                r_power_12 = np.power(r_power_6, 2)

                chain1_hydro = chain.beads[idx].hydro
                chain2_hydro = chain.beads[jdx].hydro

                if chain1_hydro=="N" or chain2_hydro=="N":
                    V_non+=4*e_h*r_power_12
                elif chain1_hydro=="B" and chain2_hydro=="B":
                    BB_lambda = np.random.normal(loc=1.0, scale=lambda_std)
                    V_non+= BB_lambda*4*e_h*(r_power_12-r_power_6)
                else:
                    V_non+=4*e_L*(r_power_12+r_power_6)

            else:
                Vnon+=100*e_h

    return V_non

def compute_energy(chain, a=1.7, e_h=1, lambda_std=10):

    Energy=0
    Energy+= compute_V_DIH(chain, e_h)
    Energy+= compute_V_BA(chain, e_h)
    Energy+= compute_V_BL(chain, a, e_h)
    Energy+= compute_V_non(chain, a, e_h, lambda_std)

    return Energy


def bead_gradient(chain, bead_idx, h_step=0.001, a=1.7, e_h=1, lambda_std=10):

    if bead_idx >= len(chain):
        return 0

    first_idx = max(0, bead_idx-3)
    last_idx = min(len(chain), bead_idx+4)

    steps = ([0.5*h_step,0,0], [0,0.5*h_step,0], [0,0,0.5*h_step])
    gradient = []

    for step in steps:
        bead_coords = chain[bead_idx].coords
        chain_plus = chain
        chain_plus[bead_idx].coords = bead_coords+step

        chain_minus = chain
        chain_minus[bead_idx].coords = bead_coords-step
        beads_plus=[]
        beads_minus=[]
        for idx in range(first_idx, last_idx):
            beads_plus.append(chain_plus[idx])
            beads_minus.append(chain_minus[idx])

        chain_plus = Chain(beads_plus)
        chain_minus = Chain(chain_minus)

        E_plus = compute_energy(chain_plus, a, e_h, lambda_std)

        E_minus = compute_energy(chain_minus, a, e_h, lambda_std)

        gradient.append((E_plus - E_minus))


    return gradient




# ============================================================================= #
#                                TEST FUNCTIONS
# ============================================================================= #

a=1.8
beads = [
    Bead(1.0, "N", a*0,0,-1),
    Bead(0.5, "N", a*1,2,10),
    Bead(1.2, "B", a*2, 0,1),
    Bead(1.0, "L", a*3,1,2),
    Bead(1.0, "B", a*4,-1,0),
    Bead(1.0, "L", a*5,-2,1),
    Bead(1.0, "B", a*6,2,2),
    Bead(1.0, "B", a*7,1,2),
    Bead(1.0, "B", a*8,0,1),
    Bead(1.0, "L", a*9,-2,-1),
    Bead(1.0, "B", a*10,0,0) ]

chain=Hydro_Chain(beads, avg_bond_length=a)

chain.fill_bond_forces()
chain.fill_angle_forces()
chain.fill_dihedral_forces()

#print(len(chain))
#print(chain.bond_forces)

#print(chain.bond_angle_forces)
print(chain.dihedral_forces)
# print(len(chain.chain_vector))
# print(len(chain.bond_angles))
# print(len(chain.dihedral_angles))


# print(bead_gradient(chain, 2))


# Implement Random Force
# Implement Momentum/Velocity in the beads
# Implement Drag coefficient
# Implement Verlet Algorithm
