from .reaction_diffu import *

class TwoComponent(RDSsystem):
    def __init__( self,F,G,D1=0,D2=0,space_size = 50,dim = 2,\
        dt=0.01,boundary = 'Neumann',init_dis= 'random',*args):
        RDSystem.__init__(
        self,reactants,
        space_size =space_size, 
        init_dis = init_dis,
        boundary = boundary, dim=dim,dt=dt)
        Reactant1 = Reactant(D1)
        Reactant2 = Reactant(D2)
        reactants = [Reactant1,Reactant2]
        #assign the reaction:
        self.F = F
        self.G = G
    
    def reaction(self):
        u = self.dis
        dudt_reac = np.zeros(np.shape(self.dis))      
        if self.dim==2: #2D        
            dudt_reac[0,:,:] +=  self.F(u[0],u[1])
            dudt_reac[1,:,:] +=  self.G(u[0],u[2])
        return dudt_reac
