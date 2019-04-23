
from .reaction_diffu import *
'''
A derived class for N-state system
A n-state system is a good physical description for thermophesis
'''

class State(Reactant):
    def __init__(self,D=0,E = 0):
        self.diffu_coe = D
        self._energy = E
    
    def info(self):
        print('Energy E=',self._energy,', Diffusion coefficien D=',self.diffu_coe);


class NStateSystem(RDSystem):
    def __init__(self,reactants=(),distribution=[],T_field = [],dt=0.01,space_size=50,boundary ='Neumann',*args):
        self.T_field = T_field
        self.dim = dim = np.size(np.shape(T_field))
        RDSystem.__init__(
            self,reactants,
            space_size = np.shape(T_field)[0], 
            boundary = boundary,dim=dim,dt=dt)
        #   and the spatical distribution   
        if distribution != []:
            self.dis = distribution  
        else:
            self.dis = np.ones(np.hstack((len(reactants),(np.shape(T_field)))))/len(reactants)
        dis_non_diffu = [np.exp(-rt._energy/self.T_field ) for rt in reactants]
        dis_non_diffu = [dis_non_diffu[i]/np.sum(dis_non_diffu,axis=0) for i in range(self._num_chem)]
        self.dis = np.array(dis_non_diffu)
        self.env = np.array(dis_non_diffu)

    def info(self):
        print('This particle have %i states:'%self._num_chem)
        for i in range(self._num_chem):
            print('*state',i,': ')
            self._chem_species[i].info()

    def reaction(self):
        T_field = self.T_field;
        u = self.dis
        dudt_reac = np.zeros(np.shape(u))
        def rate(E_i,E_f,T_field):
            rate_m = np.exp(-(E_f-E_i)/T_field)
            rate_m[np.where(rate_m>1)]=1
            return rate_m
        if self.dim==2: #2D        
            for i in range(self._num_chem):
                for j in range(self._num_chem):
                    Ei = self._chem_species[i]._energy
                    Ej = self._chem_species[j]._energy
                    dudt_reac[i,:,:] +=  rate(Ej,Ei,T_field)*u[j,:,:] - rate(Ei,Ej,T_field)*u[i,:,:]
        if self.dim==1:
            for i in range(self._num_chem):
                for j in range(self._num_chem):
                    Ei = self._chem_species[i]._energy
                    Ej = self._chem_species[j]._energy
                    dudt_reac[i,:] +=  rate(Ej,Ei,T_field)*u[j,:] - rate(Ei,Ej,T_field)*u[i,:]
        return dudt_reac

    def soret_coeff(self):
        ## soret = (<E*D>-<E><D>)/(<D>kT^2)
        T = self.T_field
        states = self._chem_species
        b_factor = lambda x: np.exp(-x/T);
        E_ave = 0
        D_ave = 0
        ave_ED = 0
        partition = 0
        num = np.size(states)
        for i in range(num):
            D,E = states[i].diffu_coe,states[i]._energy
            E_ave += E*b_factor(E)
            D_ave += D*b_factor(E)
            ave_ED += E*D*b_factor(E)
            partition += b_factor(E)
        E_ave = E_ave/partition
        D_ave = D_ave/partition
        ave_ED = ave_ED/partition
        soret_coe = (ave_ED-E_ave*D_ave)/D_ave/T**2
        return soret_coe

    def soret_compare(self):
        u = np.sum(self.dis,axis=0)
        T = self.T_field
        d = np.gradient(u)
        #print(d)
        #The numerical Soret coefficients
        if self.dim == 1:
            x = self.x
            num_resu = -np.gradient(u)/np.gradient(T)/u
            #The analytical Soret coefficients
            ana_resu = self.soret_coeff()
            #fig=plt.figure(figsize=(18, 9), dpi= 80, facecolor='w', edgecolor='k') 
            scale_font()
            plt.plot(x,num_resu,label = 'numerical')
            plt.plot(x,ana_resu,'--',label = 'analytical')
            plt.title('Soret coefficient')
            plt.xlabel('x')
            plt.ylabel('Soret coefficient')
            #plt.text(x[round(self.size/6)],(num_res[1]+num_res[-1])/2,'$T_1$ = %.1f, $T_2$ =%.1f \n $D_1$=%.1f, $D_2$=%.1f \n $\Delta E$=%.1f'%(self.T1,self.T2,self.D1,self.D2,self.dE))
            plt.legend()
            plt.grid()
        if self.dim == 2:
            x, y  = self.x, self.y
            #The analytical Soret coefficients
            ana_resu = self.soret_coeff();
            fig=plt.figure(figsize=(18,10), dpi= 80, facecolor='w', edgecolor='k')
            ax = fig.add_subplot(111, projection='3d')
            # Plot a basic wireframe.
            ax.plot_wireframe(x,y, ana_resu ,  rstride=10, cstride=10)
            num_resu = -np.gradient(u,axis=1)/np.gradient(T,axis=1)/u
            ax.plot_surface(x,y, num_resu , cmap=cm.coolwarm, linewidth=0, antialiased=False)
            plt.title('Soret Coffecient\n theory(mesh), simulation(surface)')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('S_T')

    def plot_T_field(self):
        T_field = self.T_field
        if self.dim == 1:
            x = self.x
            scale_font()
            plt.plot(x,T_field)
            plt.grid()
            plt.title('temperature field')
            plt.xlabel('x')
            plt.ylabel('temperature')
        if self.dim == 2:
            X,Y= self.x,self.y
            scale_font()
            fig=plt.figure(figsize=(18, 13), dpi= 80, facecolor='w', edgecolor='k')
            ax = fig.add_subplot(111, projection='3d')
            # Plot a basic wireframe.
            ax.plot_surface(X,Y, T_field, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            plt.title('Temperature field')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('T')

    def plot_stat_dis(self):
        #fig=plt.figure(figsize=(18, 9), dpi= 80, facecolor='w', edgecolor='k') 
        u = np.sum(self.dis,axis=0)
        if self.dim == 1:
            x = self.x
            for i in range(self._num_chem):
                linetype = (0,tuple([i,2,1,2,i*2,2]))
                pl.plot(x,self.dis[i,:],linestyle = linetype ,label = 'state'+str(i))
            pl.plot(x,u,label = 'particle')
            pl.xlabel('x')
            pl.ylabel('distribution')
            pl.legend()
            plt.grid()
            plt.title('Stationary distribution')
            scale_font()
        if self.dim == 2:
            #only polt the distribnution of particle, not states
            X,Y= self.x,self.y
            fig=plt.figure(figsize=(18, 13), dpi= 80, facecolor='w', edgecolor='k')
            ax = fig.add_subplot(111, projection='3d')
            # Plot a basic wireframe.
            ax.plot_surface(X,Y, u, cmap=cm.coolwarm, linewidth=0, antialiased=False)
            plt.title('Concentration of Particles   ')
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('C')


