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

    def __init__(self,reactants=(),space_size = 50,dim = 1,dt=0.01,init_dis= 'random',boundary = 'Neumann',*args):
        #chem-dic is a dictionary to store the...
        #   and the spatical distribution
        if len(reactants) == 0:
            print('No chemical species added')
        else:
            self.size = space_size
            self.dx = 1/space_size
            self._num_reactants = len(reactants)
            self.reactants = reactants
            self.dim = dim
            self.dt = min(1/5/space_size**2/max([reat.diffu_coe for reat in self.reactants]),dt);
            #print('dt=',self.dt)
            self.boundary = boundary
            #print(self.dt)
            if dim ==1:
                self.x = np.linspace(0,1,space_size)
                self.dis = self.generator[init_dis](self.x,dim)
            if dim==2:
                space = np.linspace(0,1,space_size)
                self.x,self.y = np.meshgrid(space,space)
                self.dis = self.generator[init_dis]([self.x,self.y],dim)
            self.dis = [self.dis /self._num_reactants for i in range(self._num_reactants)]
      
                
    def diffusion(self):
        u = self.dis
        dudt_diff = np.zeros(np.shape(u))
        for i in range(self._num_reactants):
            dudt_diff[i] = ndimage.laplace(u[i])/self.dx**2*self.reactants[i].diffu_coe;
        return dudt_diff

    def reaction(self):
        return np.zeros(np.shape(self.dis))
    
    def diffusion_reaction(self,boundary = 'Neumann'):
        self.dis += (self.diffusion()+self.reaction())*self.dt;
        if self.boundary == 'Neumann':
            for u in self.dis:
                if self.dim == 1:
                    u[0]=u[1]
                    u[-1]=u[-2]
                if self.dim ==2:
                    u[0]=u[1]
                    u[-1]=u[-2]
                    u[:,0]=u[:,1]
                    u[:,-1]=u[:,-2]
        if self.boundary == 'periodic':
            for u in self.dis:
                if self.dim == 1:
                    u[0]=u[-1]=(u[0]+u[-1])/2
                if self.dim ==2:
                    u[0]=u[-1]=(u[0]+u[-1])/2
                    u[:,0]=u[:,-1]= (u[:,0]+u[:,-1])/2
        if self.boundary == 'Dirichlet':
            for i in range(self._num_reactants):
                if self.dim == 1:
                    self.dis[i,0]=self.env[i,0]
                    self.dis[i,-1]=self.env[i,-1]
                if self.dim ==2:
                    self.dis[i,0]=self.env[i,0]
                    self.dis[i,-1]=self.env[i,-1]
                    self.dis[i,:,0]=self.env[i,:,0]
                    self.dis[i,:,1] =self.env[i,:,-1]
    '''
    def integrate(self,t=10):
        solver = ode(self.diffusion()+self.reaction()).set_integrator('vode', method='bdf', atol=1e-8, rtol=1e-8, nsteps=5000 )
        solver.set_initial_value(self.dis, 0)
        solver.integrate(t)
        self.dis = solver.y
    '''
        
    def evolution(self,evolution_time,print_time = False):
        i=0
        for i in range(int(evolution_time/self.dt)):
            self.diffusion_reaction();
            if i%10000 ==0:
                print('dt: ',i*self.dt) if print_time == True else None
            #print(np.sum(abs(self.dis-dis_tem))
        print('run time=',i*self.dt) if print_time == True else None

    def plot_dis(self):
        #fig = plt.figure()
        n = self._num_reactants
        for i in range(n):
            plt.subplot(1,n,i+1)
            plt.imshow(self.dis[i])
            plt.axis('off')
        #plt.show()

    def stationary(self,optimize_target=1e-7):
        i=0
        loss_target = optimize_target# /np.size(self.dis);
        print('target: ',loss_target)
        while(1):
            dis_tem = np.array(self.dis)
            self.diffusion_reaction();
            if i%5000==0:
                print(i,'run, loss:',np.sum(abs(self.dis-dis_tem)))
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
