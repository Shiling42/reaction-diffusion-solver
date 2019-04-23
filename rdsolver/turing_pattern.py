from .reaction_diffu import *

class TuringPattern(RDSystem):
    def __init__( self,a = 2.8e-4,b=5e-3,tau=0.1, k=-.005,space_size = 50,dim = 2,\
        dt=0.01,boundary = 'Neumann',init_dis= 'random',*args):
        self.a,self.b,self.tau,self.k = (a,b,tau,k )
        Reactant1 = Reactant(a)
        Reactant2 = Reactant(b/tau)
        reactants = [Reactant1,Reactant2]
        RDSystem.__init__(
        self,reactants,
        space_size =space_size, 
        init_dis = init_dis,
        boundary = boundary, dim=dim,dt=dt)
    
    def reaction(self):
        u = self.dis
        dudt_reac = np.zeros(np.shape(self.dis))      
        if self.dim==2: #2D        
            dudt_reac[0,:,:] +=  u[0]- u[0]**3 - u[1] + self.k
            dudt_reac[1,:,:] +=  (u[0] - u[1])/self.tau
        return dudt_reac
    
    def evolution(self,evolution_time):
        i=0
        for i in range(int(evolution_time/self.dt)):
            self.diffusion_reaction();
            if i%1000 ==0:
                print('time: ',i*self.dt)
            #print(np.sum(abs(self.dis-dis_tem))
        print('running time=',i*self.dt)


if __name__ == "__main__":   
    Turing_1 = TuringPattern()
    Turing_1.stationary()
    plt.imshow(Turing_1.dis)
    plt.show()
5