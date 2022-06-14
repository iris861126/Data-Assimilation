"""
Create the observations for DA experiment
Load:
  x_t.txt
Save:
  x_o_save.txt
"""
import numpy as np
from scipy.integrate import ode
import lorenz96
from settings import *

# load the nature run
x_t = np.genfromtxt('x_t.txt')
x_t_4d = np.genfromtxt('x_t_4d.txt')

x_o_4d_save = np.zeros((nT_4d+1,N))
x_o_save = np.zeros((nT+1,N))

for k in range(N):
    obs_err = np.random.normal(mean_obs_err,std_obs_err,nT_4d+1)
    x_o_4d_save[:,k] = x_t_4d[:,k] + obs_err
    print('mean =', np.mean(obs_err) , "\nstd varible =", np.std(obs_err))

x_o_save[0] = x_o_4d_save[0,:]
for i in range(nT_4d+1):
    if(i%4 == 0.):
        x_o_save[int(i/4.),:] = x_o_4d_save[i,:]
#%%
#test A
a = 1
for i in range(2,32):
    if(i%2 != 0):
        a = np.hstack((a,i))
x_o_density = np.delete(x_o_save,a,axis=1)
H = np.eye(N)
H = np.delete(H,a,axis=0)

#test B
a = 1
for i in range(2,32):
    if(i%3 != 0):
        a = np.hstack((a,i))
x_o_density2 = np.delete(x_o_save,a,axis=1)
H2 = np.eye(N)
H2 = np.delete(H2,a,axis=0)

#test C
b = [1,3,5,7,11,13,15,17,19,21,25,27,29,31,33,35]
x_o_density3 = np.delete(x_o_save,b,axis=1)
H3 = np.eye(N)
H3 = np.delete(H3,b,axis=0)

# save the observations for DA experiment
np.savetxt('x_o_4d_save.txt', x_o_4d_save)
np.savetxt('x_o_save.txt', x_o_save)
np.savetxt('x_o_density.txt',x_o_density)
np.savetxt('x_o_density2.txt',x_o_density2)
np.savetxt('x_o_density3.txt',x_o_density3)

np.savetxt('H.txt',H)
np.savetxt('H2.txt',H2)
np.savetxt('H3.txt',H3)



