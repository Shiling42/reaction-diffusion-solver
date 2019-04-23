
import sys
#import unittest
sys.path.append('/home/shiling/Documents/GitHub/reaction-diffusion-solver')

from rdsolver.turing_pattern import *
#from .. import rdsolver



Turing_1 = TuringPattern(space_size = 50,
        boundary = 'periodic',init_dis='random')
#Turing_1.evolution(30)
Turing_1.stationary(optimize_target=0.02)
fig = plt.figure()
plt.subplot(1,2,1)
plt.imshow(Turing_1.dis[0])
plt.subplot(1,2,2)
plt.imshow(Turing_1.dis[1])
plt.show()
