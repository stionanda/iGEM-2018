import sys, random
import pygame
from pygame.locals import *
import math
import pymunk
import pymunk.pygame_util

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation


temp1=int( raw_input('Enter the temperature (4-37 Celsius): ') ) 
temp2=temp1+273.15
Boltzmann=1.38066e-23
Protein_Mass=5e-23
Protein_Velocity= math.sqrt((3*Boltzmann*temp2)/Protein_Mass)

# particle_types
collision_types={
    "c1": 1,
    "n1": 2,
    "RFP": 3,
    "c2": 4,
    "n2": 5,
    "GFP": 6,    
}

RED = (255,0,0)
GREEN = (0,255,0)
BLUE = (0,0,255)
BLACK= (0,0,0)
WHITE = (255,255,255)
MAGENTA = (255,0,255) 
CYAN = (0,255,255)
YELLOW = (255,255,0)
ORANGE = (255,165,0)

R=4
R1=10
D=10.0
S=500

class World:
    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.molecules = []
        self.particles = []
        self.space = pymunk.Space()
        self.space.gravity = (0.0, 0.0)
        self.add_border()
        self.joints =[]
        return 

    def add_border(self):
        body = pymunk.Body(body_type = pymunk.Body.STATIC) 
        body.position = (0, 0)
        l1 = pymunk.Segment(body, (0, 0), (self.width, 0), 5)
        l2 = pymunk.Segment(body, (0, 0), (0, self.height), 5) 
        l3 = pymunk.Segment(body, (0, self.height), (self.width, self.height), 5) 
        l4 = pymunk.Segment(body, (self.width, self.height), (self.width, 0), 5) 
        for line in l1,l2,l3,l4:
            line.elasticity  = 1.0
            line.friction = 0.0
        self.space.add(l1, l2, l3, l4) 
        return    

    def add_particle(self):
        particle = Particle(self)
        particle1 = Particle1(self)
        self.molecules.append(particle)
        self.molecules.append(particle1)
        return    


    def shape_to_particle(self,shape):
        for particle in self.particles:
            if particle.shape == shape:
                return particle
        for particle1 in self.particles:
            if particle1.shape == shape:
                return particle1
        return None
    # def checkParticle(self):

    def particle_to_molecule(self, in_particle):
        for molecule in self.molecules:
            for particle in molecule.particles:
                if particle is in_particle:
                    return molecule
        return None

class Bond:
    def __init__(self, a, b):
        self.LHS=a
        self.RHS=b
        d1 = pymunk.DampedSpring(a.body, b.body, (R,R),(-R,R), 0, S, D)
        d2 = pymunk.DampedSpring(a.body, b.body, (R,R/2),(-R,R/2), 0, S, D)
        d3 = pymunk.DampedSpring(a.body, b.body, (R,-R/2),(-R,-R/2), 0, S, D)
        d4 = pymunk.DampedSpring(a.body, b.body, (R,-R),(-R,-R), 0  , S,  D)        
        #d1.collide_bodies=False 
        self.joints=[d1,d2,d3,d4]
        self.LHS.right_bond=self
        self.RHS.left_bond=self
        a.world.space.add(d1,d2,d3,d4)

    def unbond(self):
        for joint in self.joints:
            self.RHS.world.space.remove(joint)
            joint.collide_bodies=True 
        self.LHS.right_bond=None
        self.RHS.left_bond=None

class Molecule:
    def __init__(self, world, type="RFP"):
        self.particles = []
        self.world = world
        self.world.molecules.append(self)
        if type=="RFP":
            self.make_RFP()
        if type=="GFP":
            self.make_GFP()
        return

    def make_RFP(self):
        n1=Particle(self.world,species="n1")
        RFP=Particle1(self.world,species="RFP")
        c1=Particle(self.world,species="c1")

        x0,y0 = RFP.body.position
        RFP.setPosition(0,0)
        n1.setPosition(2*R1,0)
        c1.setPosition(-2*R1,0)
       
        RFP.offsetPosition(x0,y0)
        n1.offsetPosition(x0,y0)
        c1.offsetPosition(x0,y0)
        RFP.shape.color=RED
        n1.shape.color=CYAN
        c1.shape.color=YELLOW
    
        RFP.shape.collision_type=collision_types["RFP"]
        n1.shape.collision_type=collision_types["n1"]
        c1.shape.collision_type=collision_types["c1"]

        random_impulse = (random.randint(-int(Protein_Velocity),int(Protein_Velocity)),random.randint(-int(Protein_Velocity),int(Protein_Velocity)))
        RFP.body.apply_impulse_at_local_point(random_impulse)
        n1.body.apply_impulse_at_local_point(random_impulse)
        c1.body.apply_impulse_at_local_point(random_impulse)
        
        random_spin = random.randint(-150,150)
        n1.body.apply_impulse_at_local_point((0,random_spin))
        c1.body.apply_impulse_at_local_point((0,-random_spin))

        Bond(c1,RFP)
        Bond(RFP,n1)

        self.particles.append(c1)
        self.particles.append(RFP)
        self.particles.append(n1)
        return
    
    def make_GFP(self):
        n2=Particle(self.world,species="n2")
        GFP=Particle1(self.world,species="GFP")
        c2=Particle(self.world,species="c2")

        x0,y0 = GFP.body.position
        GFP.setPosition(0,0)
        n2.setPosition(2*R1,0)
        c2.setPosition(-2*R1,0)
       
        GFP.offsetPosition(x0,y0)
        n2.offsetPosition(x0,y0)
        c2.offsetPosition(x0,y0)
        GFP.shape.color=GREEN
        n2.shape.color=ORANGE
        c2.shape.color=BLUE
    
        GFP.shape.collision_type=collision_types["GFP"]
        n2.shape.collision_type=collision_types["n2"]
        c2.shape.collision_type=collision_types["c2"]

        random_impulse = (random.randint(-int(Protein_Velocity),int(Protein_Velocity)),random.randint(-int(Protein_Velocity),int(Protein_Velocity)))
        GFP.body.apply_impulse_at_local_point(random_impulse)
        n2.body.apply_impulse_at_local_point(random_impulse)
        c2.body.apply_impulse_at_local_point(random_impulse)
        
        random_spin = random.randint(-150,150)
        n2.body.apply_impulse_at_local_point((0,random_spin))
        c2.body.apply_impulse_at_local_point((0,-random_spin))

        Bond(c2,GFP)
        Bond(GFP,n2)

        self.particles.append(c2)
        self.particles.append(GFP)
        self.particles.append(n2)
        return

class Particle:
    def __init__(self, world, mass=.5,species=None):
        self.species=species
        self.world = world
        self.left_bond = None
        self.right_bond = None
        self.mass = .5
        moment = pymunk.moment_for_box(self.mass,(2*R,2*R)) 
        self.body = pymunk.Body(self.mass, moment) 
        self.shape = pymunk.Poly(self.body, [(-R,-R),(-R,R),(R,R),(R,-R)])

        self.shape.elasticity  = 1.0
        self.shape.friction = 0.0
        self.world.space.add(self.body,self.shape)
        self.setRandomPosition()
        self.world.particles.append(self)


    def setRandomPosition(self):
        border = 3*R
        max_tries=10
        while max_tries>0:
            x=random.randint(border,self.world.width-border)
            y=random.randint(border,self.world.height-border)
            filter = pymunk.ShapeFilter(mask=pymunk.ShapeFilter.ALL_MASKS)
            results = self.world.space.bb_query(
                pymunk.BB(x-border,y-border,x+border,y+border),
                        filter)
            if len(results)==0:
                max_tries=0
        self.body.position = (x,y)
        return

    def setPosition(self,x,y):
        self.body.position = (x,y)

    def offsetPosition(self,x,y):
        x0,y0=self.body.position
        self.body.position=(x+x0,y+y0)

class Particle1:
    def __init__(self, world, mass=1.5,species=None):
        self.species=species
        self.world = world
        self.left_bond = None
        self.right_bond = None
        self.mass = 1.5
        moment = pymunk.moment_for_box(self.mass,(2*R1,2*R1)) 
        self.body = pymunk.Body(self.mass, moment) 
        self.shape = pymunk.Poly(self.body, [(-R1,-R1),(-R1,R1),(R1,R1),(R1,-R1)])
        #self.shape = pymunk.Circle(self.body, R1)

        self.shape.elasticity  = 1.0
        self.shape.friction = 0.0
        self.world.space.add(self.body,self.shape)
        self.setRandomPosition()
        self.world.particles.append(self)


    def setRandomPosition(self):
        border = 3*R
        max_tries=10
        while max_tries>0:
            x=random.randint(border,self.world.width-border)
            y=random.randint(border,self.world.height-border)
            filter = pymunk.ShapeFilter(mask=pymunk.ShapeFilter.ALL_MASKS)
            results = self.world.space.bb_query(
                pymunk.BB(x-border,y-border,x+border,y+border),
                        filter)
            if len(results)==0:
                max_tries=0
        self.body.position = (x,y)
        return

    def setPosition(self,x,y):
        self.body.position = (x,y)

    def offsetPosition(self,x,y):
        x0,y0=self.body.position
        self.body.position=(x+x0,y+y0)

def main():
    pygame.init()
    xmax=1200
    ymax=800
    screen = pygame.display.set_mode((xmax, ymax))
    clock = pygame.time.Clock()
    world = World(xmax, ymax)

    def collision_begin_RFP_n1_c2_GFP(arbiter, space, data):
        n1 = world.shape_to_particle(arbiter.shapes[0]) 
        c2 = world.shape_to_particle(arbiter.shapes[1])
        if n1.left_bond==None or c2.right_bond == None: 
            return True
        RFP_n1 =  n1.left_bond
        c2_GFP = c2.right_bond
        RFP=RFP_n1.LHS
        GFP=c2_GFP.RHS
        # 
        RFP_n1.unbond()
        c2_GFP.unbond()

        Bond(n1,c2)
        Bond(RFP,GFP)
        
        n1_mol = world.particle_to_molecule(n1)
        c2_mol = world.particle_to_molecule(c2)        
        n1_mol.particles=n1_mol.particles[:-1]+c2_mol.particles[1:]
        c2_mol.particles=[n1,c2]
        return True
    
    def collision_begin_GFP_n2_c1_RFP(arbiter, space, data):
        n2 = world.shape_to_particle(arbiter.shapes[0]) 
        c1 = world.shape_to_particle(arbiter.shapes[1])

        if n2.left_bond==None or c1.right_bond == None: 
            return True
        GFP_n2 =  n2.left_bond
        c1_RFP = c1.right_bond
        GFP=GFP_n2.LHS
        RFP=c1_RFP.RHS
        # 
        GFP_n2.unbond()
        c1_RFP.unbond()

        Bond(n2,c1)
        Bond(GFP,RFP)
        n2_mol = world.particle_to_molecule(n2)
        c1_mol = world.particle_to_molecule(c1)        
        n2_mol.particles=n2_mol.particles[:-1]+c1_mol.particles[1:]
        c1_mol.particles=[n2,c1]
        return True
    

    h1 = world.space.add_collision_handler(
        collision_types["n1"],
        collision_types["c2"])
    h1.begin = collision_begin_RFP_n1_c2_GFP
    
    h2 = world.space.add_collision_handler(
        collision_types["n2"],
        collision_types["c1"])
    h2.begin = collision_begin_GFP_n2_c1_RFP


    draw_options = pymunk.pygame_util.DrawOptions(screen)
    draw_options.flags = pymunk.SpaceDebugDrawOptions.DRAW_SHAPES

    nRFP = 100
    for i in range(0,nRFP):
        new_particle = Molecule(world,"RFP")
        
    nGFP = 100
    for i in range(0,nGFP):
        new_particle = Molecule(world,"GFP")
    mytime=0
    f=open("outdata.csv","w")
    while True:
        for event in pygame.event.get():
            if event.type == QUIT:
                f.close()
                sys.exit(0)
                break
            elif event.type == KEYDOWN and event.key == K_ESCAPE:
                f.close()
                sys.exit(0)
                break

        world.space.step(1/50.0)
        screen.fill((255,255,255))
        world.space.debug_draw(draw_options)
        pygame.display.flip()
        clock.tick(50)
        mytime+=1
        if mytime%10==0:
            n=0
            for mol in world.molecules:
                if len(mol.particles)>2:
                    n+=1
            print mytime,n
            f.write("{},{}\n".format(mytime,n))
        if mytime%10000==0:
            f.close()
            sys.exit(0)
            break
if __name__ == '__main__':
    main()
    

