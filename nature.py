"""
Generate the nature run
Save:
  x_t.txt
"""
import numpy as np
from scipy.integrate import ode
import lorenz96
from settings import *

# settings of spin-up
sigma_x0 = 0.1  # size of initial perturpation
Tspinup = 100.  # initial spin-up time

# spin-up from a random initail value
x_t_0 = sigma_x0 * np.random.randn(N)

solver = ode(lorenz96.f).set_integrator('dopri5', nsteps=10000)
solver.set_initial_value(x_t_0, 0.).set_f_params(F)
solver.integrate(Tspinup)
x_t_4d_save = np.array([solver.y], dtype='f8')

# create nature
solver = ode(lorenz96.f).set_integrator('dopri5')
solver.set_initial_value(x_t_4d_save[0], 0.).set_f_params(F)

tt_4d = 1
while solver.successful() and tt_4d <= nT_4d:
    solver.integrate(solver.t + dT_4d)
    x_t_4d_save = np.vstack([x_t_4d_save, [solver.y]])
#    print('timestep =', tt, round(solver.t, 10))
    tt_4d += 1
    
x_t_save = np.zeros((nT+1,N))
x_t_save[0] = x_t_4d_save[0,:]
for i in range(nT_4d+1):
    if(i%4 == 0.):
        x_t_save[int(i/4.),:] = x_t_4d_save[i,:]
  
# save data
np.savetxt('x_t_4d.txt', x_t_4d_save)
np.savetxt('x_t.txt', x_t_save)
