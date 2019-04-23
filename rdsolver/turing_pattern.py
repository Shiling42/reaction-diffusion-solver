from .reaction_diffu import *

class TuringPattern(RDSystem):
    def __init__(self,a = 2.8e-4,b=5e-3,tau=0.1,k=-.005,space_size = 50,dim = 2,dt=0.01,init_dis= 'random',*args):
        self.a,self.b,self.tau,self.k = (a,b,tau,k )
        Reactant1 = Reactant(a)
        Reactant2 = Reactant(b/tau)
        reactants = [Reactant1,Reactant2]
        RDSystem.__init__(
        self,reactants,
        space_size =space_size, 
        dim=dim,dt=dt)
        print(1)
    def reaction(self):
        u = self.dis
        dudt_reac = np.zeros(np.shape(u))      
        if self.dim==2: #2D        
            dudt_reac[0,:,:] +=  u[0]- u[0]**3 - u[1] + self.k
            dudt_reac[1,:,:] +=  (u[0] - u[1])*self.size**2
        return dudt_reac
    def evolution(self):
        i=0
        for i in range(int(20/self.dt)):
            dis_tem = np.array(self.dis)
            self.diffusion_reaction();
            i += 1
            #print(np.sum(abs(self.dis-dis_tem)))
            if np.sum(abs(self.dis-dis_tem)) <1e-6:
                break;
        print('running time=',i*self.dt)


if __name__ == "__main__":   
    Turing_1 = TuringPattern()
    Turing_1.stationary()
    plt.imshow(Turing_1.dis)
    plt.show()
