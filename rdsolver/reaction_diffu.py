import pylab as pl
import numpy as np
from numpy import linalg as LA
# impot modules from sicpy
from scipy import ndimage, misc
from scipy.linalg import block_diag
import scipy.constants as sc
import scipy as sp
from scipy.integrate import ode
# a module for periodically choose line style
from itertools import cycle
# plot mopdules
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
        #dudt_diff = np.zeros(np.shape(u))
        if self.boundary == 'Dirichlet':
            dudt_diff = [ndimage.laplace(u[i]*self.reactants[i].diffu_coe,mode='constant')/self.dx**2 for i in range(self._num_reactants)]
        if self.boundary == 'Neumann':
            dudt_diff = [ndimage.laplace(u[i]*self.reactants[i].diffu_coe,mode='nearest')/self.dx**2 for i in range(self._num_reactants)]
        if self.boundary == 'periodic':
            dudt_diff = [ndimage.laplace(u[i]*self.reactants[i].diffu_coe,mode='wrap')/self.dx**2 for i in range(self._num_reactants)]
        return dudt_diff

    def reaction(self):
        return np.zeros(np.shape(self.dis))

    def noise(self):
        return np.zeros(np.shape(self.dis))
    
    def diffusion_reaction(self):
        self.dis += (self.diffusion()+self.reaction()+self.noise()/np.sqrt(self.dt))*self.dt;

      
    def evolve(self,evolve_time,print_time = False,draw=False):
        i=0
        n_r = self._num_reactants
        if draw == True:
            fig, axs = plt.subplots(1, n_r, figsize=(7, 8*n_r), sharey=True)
        for i in range(int(evolve_time/self.dt)):
            self.diffusion_reaction();
            if i%500 ==0:
                print('dt: ',i*self.dt) if print_time == True else None
                if draw == True:
                    for i in range(n_r):
                        axs[i].imshow(self.dis[i])
                        axs[i].set_title('Reactant %i'%i)
                        axs[i].axis('off')
                    plt.draw()
                    plt.pause(0.00000001)
            #print(np.sum(abs(self.dis-dis_tem))
        print('run time=',i*self.dt) if print_time == True else None

    def plot_dis(self):
        n = self._num_reactants
        fig, axs = plt.subplots(1, n, figsize=(7, 8*n), sharey=True)
        for i in range(n):
            axs[i].imshow(self.dis[i])
            axs[i].set_title('Reactant %i'%i)
            axs[i].axis('off')
        #fig.suptitle('Concentration distribution')
        plt.plot()

    def stationary(self,optimize_target=1e-6,print_time = False):
        i=0
        loss_target = optimize_target# /np.size(self.dis);
        print('target: ',loss_target) if print_time == True else None
        while(1):
            dis_tem = np.array(self.dis)
            self.diffusion_reaction();
            if print_time:
                if i%5000==0 :
                    print(i,'run, loss:',np.sum(abs(self.dis-dis_tem)))
                    print('target: ',loss_target)
            i += 1
            #print(np.sum(abs(self.dis-dis_tem)))
            if np.sum(abs(self.dis-dis_tem)) <loss_target:
                break;
        print('running time=',i*self.dt) if print_time == True else None
        #return self.dis


