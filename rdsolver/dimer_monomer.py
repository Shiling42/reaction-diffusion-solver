from .n_state_system import *

class MonomerDimer(NStateSystem):
    def __init__(self, D_mon = 2, D_dimer = 1,
                ,space_size = 30,dim = 2, dt=0.01,distribution=[],T_field = [],boundary = 'Neumann',init_dis= 'random',*args):
        self.a,self.b,self.tau,self.k = (a,b,tau,k )
        reactant1 = Reactant(D_dimer)
        reactant2 = Reactant(D_monomer)
        reactants = [Reactant1,Reactant2]
        NStateSystem.__init__(self,reactants=(reactant1,reactant2),distribution=[],T_field = T_field,boundary ='Neumann',*args)
        """
        RDSystem.__init__(
        self,reactants,
        space_size =space_size, 
        init_dis = init_dis,
        boundary = boundary, dim=dim,dt=dt)
        """
    
    def reaction(self):
        T_field = self.T_field;
        u = self.dis
        dudt_reac = np.zeros(np.shape(u))
        k12 = k0 * np.maximum(np.exp(-(E_monomer-E_dimer)/T_field),1)
        k21 = k0 * np.maximum(np.exp(-(E_dimer-E_monomer)/T_field),1)
        dudt_reac[0] += + 2 * (k12 * u[1] - k21 * np.power[u[0],2]);
        dudt_reac[1] += - 1/2 * dudt_reac[0] 
        return dudt_reac

if __name__ == "__main__":   
    Turing_1 = TuringPattern()
    #Turing_1.stationary()
    t=np.linspace(0,10,1000)
    Turing_1.integrate(t)
    plt.imshow(Turing_1.dis[0])
    plt.show()

