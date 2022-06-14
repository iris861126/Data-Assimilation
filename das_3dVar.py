"""
The data assimilation system: incremental 3dVar
Load:
  x_a_init.txt
Save:
  x_b_in3d.txt
  x_a_in3d.txt
"""
import numpy as np
from scipy.integrate import ode
import lorenz96
from settings import *
import math
import scipy.optimize as op

#the cost func.
def costFunction(x_g):
    y_b = np.dot(H, x_b)
    d = y_o - y_b
    dx = x_g-x_b
    J_3d = 0.5*np.dot(np.dot(np.transpose(dx),np.linalg.inv(Pb_3dVar)),(dx)) + \
           0.5*np.dot(np.dot(np.transpose(np.dot(H,dx)-d),np.linalg.inv(R)),(np.dot(H,dx)-d))
    return(J_3d)

#load initial condition
x_a_init = np.genfromtxt('x_a_init.txt')

# load observations
y_o_save = np.genfromtxt('x_o_save.txt')

# load truth
x_t_save = np.genfromtxt('x_t.txt')

# initial x_b: no values at the initial time (assign NaN)
x_b_save = np.full((1,N), np.nan, dtype='f8')

# initial x_a: from x_a_init
x_a_save = np.array([x_a_init])

#load initial BEC
Pb_3dVar = np.genfromtxt('Pb.txt')    #initial BEC

# load H
#H = np.genfromtxt('H2.txt')
H = np.eye(N)

# calculate for obs. error covariance(R)

R = np.zeros((N,N))
for i in range(N):
     R[i,i] = (std_obs_err)**2.
"""
# density obs. err 
R = np.zeros((19,19))
for i in range(19):
     R[i,i] = (std_obs_err)**2.
"""
#===========================================
#DA process: incremental 3dVar
tt = 1    
while tt <= nT:
    print('3dVar DA process...')
    tts = tt - 1
    Ts = tts * dT  # forecast start time
    Ta = tt  * dT  # forecast end time (DA analysis time)
    print('Cycle =', tt, ', Ts =', round(Ts, 10), ', Ta =', round(Ta, 10))
    
    #--------------
    # forecast step
    #--------------
    
    solver = ode(lorenz96.f).set_integrator('dopri5')
    solver.set_initial_value(x_a_save[tts], Ts).set_f_params(F)
    solver.integrate(Ta)
    x_b_save = np.vstack([x_b_save, [solver.y]])
    
    #--------------
    # analysis step
    #--------------
        
    #truth
    x_t = x_t_save[tt].transpose()
    # background
    x_b = x_b_save[tt].transpose()
    
    # observation
    y_o = y_o_save[tt].transpose()
    
    iter = 0
    x_g = x_b    #first guess for xa[x_g]
    x_g_old = x_g
    while True:
        #minimize the cost func.
        result = op.minimize(fun=costFunction,x0=x_g,method='CG')
        x_g = result.x
        J_3d = result.fun
        if (iter < 5):         
            x_g_old = x_g
            #print(x_g)
            print(J_3d)
            iter += 1
        else:
            x_a = x_g
            print('x_a= ',x_a,'Minimizing done!!')
            break
    
    
    x_a_save = np.vstack([x_a_save, x_a.transpose()])
    tt += 1  


# save incremental 3dVar background and analysis data

np.savetxt('x_b_in3d.txt', x_b_save)
np.savetxt('x_a_in3d.txt', x_a_save)
"""
np.savetxt('x_b_in3d_obs1.txt', x_b_save)
np.savetxt('x_a_in3d_obs1.txt', x_a_save)
"""