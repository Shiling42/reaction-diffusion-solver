
import sys
sys.path.append('/home/shiling/Documents/GitHub/reaction-diffusion-solver')


from rdsolver import turing_pattern as tp
#dir(rds)
from rdsolver.n_state_system import *

# Define the states
D1 = 8e-4
D2 = 1.1*D1
state1 = State(D=D1,E=2)
state2 = State(D=D2,E=3.4)
states = (state1,state2)

size = 30
x = np.linspace(0,1,size)
T1,T2 = 1.6, 1
T_field = 2*T1*(1-x)*x+T2
T_field = T1+(T2-T1)*x
system1=NStateSystem(states,T_field=T_field,space_size=size)
system1.stationary()
fig=plt.figure(figsize=(18, 13), dpi= 80, facecolor='w', edgecolor='k')
plt.subplot(3,1,1)
system1.plot_T_field()
plt.subplot(3,1,2)
system1.plot_stat_dis()
plt.subplot(3,1,3)
system1.soret_compare()
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace= 0.7)
plt.show()
'''
# creat a particle
state1 = State(D=1e-3,E=2)
state2 = State(D=2e-3,E=3.4)
states = (state1,state2)

#assign a Temperature field
size = 30
x = np.linspace(-0.5,0.5, size)
y = np.linspace(-0.5,0.5, size)
xv, yv = np.meshgrid(x, y)
T1=1.6
T_field = np.exp(-(xv**2+yv**2)/0.4)*np.sqrt(1/np.pi/0.4)*T1

#create an instance of class
system2=NStateSystem(states,T_field=T_field)

system2.stationary()

# visualize
system2.plot_T_field()

system2.soret_compare()

system2.plot_stat_dis()
plt.show()
## plot the tempreature field
#system2.plot_T_field()
'''
