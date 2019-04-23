import pylab as pl
import numpy as np
from numpy import linalg as LA

from scipy import ndimage, misc
from scipy.linalg import block_diag
import scipy.constants as sc
import scipy as sp
from scipy.integrate import ode


from itertools import cycle

import random

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

from .plot_style import *



class Reactant:
    def __init__(self,diff_coe = 0):
        self.diffu_coe= diff_coe


class RDSystem:
    generator = {
            'gaussian' : (lambda x,dim, sigma=0.4: 
                np.exp(-[sum((v-np.average(v))**2 for v in x) 
                if dim==2 else (x-np.average(x))**2][0]
                /(2*sigma**2))*np.sqrt(1/(np.pi*sigma*2))),
            'uniform' : (
                lambda x,dim:
                np.ones(np.shape([sum(v for v in x) if dim==2 else x][0]))),
            'random' : (
                lambda x,dim: 
                np.random.random_sample(
                    np.shape([sum(v for v in x) if dim==2 else x][0])))
        }
    def __init__(self,reactants=(),space_size = 50,dim = 1,dt=0.01,init_dis= 'random',boundary = 'Newmann',*args):
        #chem-dic is a dictionary to store the...
        #   and the spatical distribution
        if len(reactants) == 0:
            print('No chemical species added')
        else:
            self.size = space_size
            self._num_chem = len(reactants)
            self._chem_species = reactants
            self.dim = dim
            self.dt = min(1/3/space_size**2/np.max([i.diffu_coe for i in reactants]),dt);
            print('dt=',self.dt)
            self.boundary = boundary
            #print(self.dt)
            if dim ==1:
                self.x = np.linspace(0,1,space_size)
                self.dis = self.generator[init_dis](self.x,dim)
            if dim==2:
                space = np.linspace(0,1,space_size)
                self.x,self.y = np.meshgrid(space,space)
                self.dis = self.generator[init_dis]([self.x,self.y],dim)
            self.dis = [self.dis for i in range(self._num_chem)]

    @staticmethod        
    def Laplacian(u,dim,dx,boundary='Newmann'):
        if boundary == 'Newmann':
            if dim ==1:#Firstly for 1D
                uright = u[2:]
                uleft = u[:-2]
                ucenter = u[1:-1]
                return np.hstack((u[1]-u[0],(uleft  + uright -2 * ucenter),u[-2]-u[-1]))/dx**2
            if dim == 2:            # then for 2D system
                shape = np.shape(u)
                inboundary = lambda i,j: (i>=0)&(i<shape[0])&(j>=0)&(j<shape[1])
                d2ud2x = np.zeros(np.shape(u))/dx**2
                for i in range(shape[0]):
                    for j in range(shape[1]): 

                        d2ud2x[i][j]=(
                            ((u[i][j+1]-u[i][j]) if inboundary(i,j+1)==1 else 0)
                            +((u[i][j-1]-u[i][j]) if inboundary(i,j-1)==1 else 0)
                            +((u[i-1][j]-u[i][j]) if inboundary(i-1,j)==1 else 0)
                            +((u[i+1][j]-u[i][j]) if inboundary(i+1,j)==1 else 0)
                            )/dx**2
                return d2ud2x
        if boundary == 'periodic':
            if dim ==1:#Firstly for 1D
                uright = np.roll(u,1)
                uleft = np.roll(u,-1)
                ucenter = u
                return (uright  + uleft - 2 * ucenter)/dx**2
            if dim == 2:            # then for 2D system
                utop = np.roll(u,-1,axis=0)
                ubottom = np.roll(u,1,axis=0)
                uleft = np.roll(u,1,axis=1)
                uright = np.roll(u,-1,axis=1)
                ucenter = u
                return (utop + uleft + ubottom + uright - 4 * ucenter)/dx**2
                
    def diffusion(self):
        u = self.dis
        chem_species = self._chem_species
        dudt_diff = np.zeros(np.shape(u))
        ## 1D
        if self.dim==1:
            for i in range(self._num_chem ):
                dudt_diff[i] = self.Laplacian(u[i],self.dim,1/self.size,self.boundary)*chem_species[i].diffu_coe;
        ## 2D
        if self.dim==2:
            for i in range(self._num_chem ):
                dudt_diff[i] = self.Laplacian(u[i],self.dim,1/self.size,self.boundary)*chem_species[i].diffu_coe;        
        return dudt_diff

    def reaction(self):
        u = self.dis
        return np.zeros(np.shape(u))
    
    def diffusion_reaction(self,boundary = 'Newmann'):
        self.dis += (self.diffusion()+self.reaction())*self.dt;
        #if self.boundary == 'open'
        #self.dis
    
    def integrate(self,t=1):
        solver = ode(self.diffusion_reaction()).set_integrator('vode', method='bdf', atol=1e-8, rtol=1e-8, nsteps=5000 )
        solver.set_initial_value(self.dis, 0)
        solver.integrate(t)
        self.dis = solver.y

    def stationary(self):
        i=0
        loss_target = 1e-7
        print('target: ',loss_target)
        while(1):
            dis_tem = np.array(self.dis)
            self.diffusion_reaction();
            if i%1000==0:
                print(i,'run, loos:',np.sum(abs(self.dis-dis_tem)))
                print('target: ',loss_target)
            i += 1
            #print(np.sum(abs(self.dis-dis_tem)))
            if np.sum(abs(self.dis-dis_tem)) <loss_target:
                break;
        print('running time=',i*self.dt)
        #return self.dis
    

'''
if __name__ == "__main__":   
    
    Turing_1 = TuringPattern(
        a=2.8e-4,b=5e-3,tau=.1,k=-.005,
        space_size=50,dt=0.0005)
    Turing_1.evolution()
    dist =Turing_1.dis[0]
    plt.imshow(Turing_1.dis[0])
    plt.show()
'''

if __name__ == "__main__":   
    state1 = State(D=1,E=2)
    state2 = State(D=1.3,E=3.4)
    states = (state1,state2)
    #assign Temperature field
    size = 30
    x = np.linspace(-0.5,0.5, size)
    y = np.linspace(-0.5,0.5, size)
    xv, yv = np.meshgrid(x, y)
    T1=1.6
    sigma = 0.1
    T_field = np.exp(-(xv**2+yv**2)/sigma)*np.sqrt(1/np.pi/sigma)*T1
    #create an instance of class
    system2=NStateSystem(states,T_field=T_field)
    ## plot the tempreature field
    #system2.plot_T_field()
    system2.stationary()
    system2.soret_compare()
    plt.show()
    system2.plot_T_field()
    plt.show()
    system2.plot_stat_dis()
    plt.show()
