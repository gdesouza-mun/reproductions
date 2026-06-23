from chains_and_beads import *

#def to impose periodic boundaries on a bead

class Box():
    def __init__(self, size=60, fric_coeff=1, kbT=4.2e-21):
        self.size=size
        self.fric_coeff=fric_coeff
        self.center=[size/2, size/2, size/2]
        self.kbT=kbT

        self.rng=np.random.default_rng()


    def periodic_boundaries_on_bead(self, bead):
        bead.coords = bead.coords % self.size

    def periodic_boundaries_on_chain(self, chain):
        for bead in chain:
            self.periodic_boundaries_on_bead(bead)

    def centralize_chain(self, chain):
        chain.centralize(box.center)

    def get_drag(self, chain):
        drag_forces=[]
        for bead in chain.beads:
            drag_forces.append(self.fric_coeff*bead.velocity)

        return drag_forces
    def get_noise(self, chain, step_size=1e-11):

        std = np.sqrt(2*self.fric_coeff*self.kbT/step_size)

        white_noise = self.rng.normal(0.0, scale=std, size=(len(chain),3))

        return white_noise


#Integrate

box = Box()
#Def to impose periodi cbondaries on a chain
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

rng = np.random.default_rng()
for bead in chain:
    bead.velocity = 2*(rng.random(3))-1


box.centralize_chain(chain)


print(box.get_noise(chain))
