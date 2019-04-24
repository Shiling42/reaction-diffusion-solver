
import sys
#import unittest
sys.path.append('/home/shiling/Documents/GitHub/reaction-diffusion-solver')

from rdsolver.turing_pattern import *
#from .. import rdsolver

Turing_1 = TuringPattern(space_size = 50,
        boundary = 'periodic',init_dis='random')
fig = plt.figure()
print('t=',0)
Turing_1.plot_dis()
for i in range(3):
    delta_t = 5
    Turing_1.evolution(delta_t)
    #Turing_1.stationary(optimize_target=0.05)
    print('t=',delta_t*(i+1)
    fig = plt.figure()
    Turing_1.plot_dis()
    #plt.title('t=%.1f'%t)
    #plt.axis('off')
plt.show()
