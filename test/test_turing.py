
#import sys
#import unittest
#sys.path.append('/home/shiling/Documents/GitHub/reaction-diffusion-solver')

import sys
sys.path.append('..')

from rdsolver.turing_pattern import *
#from .. import rdsolver


Turing_1 = TuringPattern(space_size = 60,
        boundary = 'periodic',init_dis='random')
#fig = plt.figure()
print('t=',0)

delta_t = 100
Turing_1.evolve(delta_t,print_time = True,draw=True)


