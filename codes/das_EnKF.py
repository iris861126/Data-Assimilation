"""
The data assimilation system:EnKF
Load:
  x_a_init.txt
Save:
  x_b_EnKF.txt
  x_a_EnKF.txt
"""
import numpy as np
from scipy.integrate import ode
import lorenz96
from settings import *
import math

def Dis(x1,x2,Lc):
    if(abs(x1-x2)>20):
       return (N-abs(x1-x2))
    else:
       return (abs(x1-x2))

# load initial condition
x_a_init = np.genfromtxt('x_a_init.txt')
x_a_en_initial = np.zeros((N,en_member))
en_initial_err = np.zeros((N,en_member))

#generate the initial condi. for ensemble mambers
for i in range(N):  #40個x
    en_initial_err[i,:] = np.random.normal(mean_en_initial_err,std_en_initial_err[i],en_member)    
for j in range(en_member):   #每個x的100個members
    x_a_en_initial[:,j] = x_a_init + en_initial_err[:,j]

np.savetxt('x_a_EnKF_initial.txt', x_a_en_initial)
# load observations
y_o_save = np.genfromtxt('x_o_save.txt')

# initial x_b: no values at the initial time (assign NaN)
x_b_en_save = np.zeros((nT+1,N,en_member))
x_b_en_save[0] = np.nan 

# initial x_a: from x_a_init
x_a_en_save = np.full((1,N,en_member), np.nan, dtype='f8')
x_a_en_save[0] = x_a_en_initial 

# load truth
x_t_save = np.genfromtxt('x_t.txt')

# load H
#H = np.genfromtxt('H2.txt')
H = np.eye(N)

# calculate for obs. error covariance(R)

R = np.zeros((N,N))
for i in range(N):
     R[i,i] = (std_obs_err)**2.
"""
R = np.zeros((19,19))
for i in range(19):
     R[i,i] = (std_obs_err)**2.
"""
#===========================================
# DA process: EnKF
tt = 1    
while tt <= nT:
    print('EnKF DA process...')
    tts = tt - 1
    Ts = tts * dT  # forecast start time
    Ta = tt  * dT  # forecast end time (DA analysis time)
    print('Cycle =', tt, ', Ts =', round(Ts, 10), ', Ta =', round(Ta, 10))
    
    #--------------
    # forecast step
    #--------------
    x_b_each = np.zeros((N,en_member)) 
    for k in range(en_member):
        solver = ode(lorenz96.f).set_integrator('dopri5')
        solver.set_initial_value(x_a_en_save[tts,:,k], Ts).set_f_params(F)
        solver.integrate(Ta)
        x_b_each[:,k] = solver.y
    x_b_en_save[tt] = x_b_each

    #-------------- 
    # analysis step
    #--------------
    
    # background
    x_b_mean = np.mean(x_b_each,axis=1)
    
    # observation
    y_o = y_o_save[tt].transpose()
    y_o_perturb = np.zeros((N,en_member))   
    y_o_en = np.zeros((N,en_member))
    for i in range(N):  #其中一個y
        y_o_perturb[i,:] = np.random.normal(0.,std_obs_err,en_member)    
    for j in range(en_member):   #這個y的100個members
        y_o_en[:,j] = y_o + y_o_perturb[:,j]

    # innovation
    y_b = H.dot(x_b_each)
    y_b_mean = np.mean(y_b,axis=1)
    
    
    #B-localization
    c_bloc = np.zeros((N,N))
    for l in range(N):
        for k in range(N): 
            if(abs(k-l)>20):
                dis = N-abs(k-l)
            else:
                dis = abs(k-l)
            if(dis>3.65*Lc):
                c_bloc[l,k] = 0.
            else:
                c_bloc[l,k] = np.exp(-(dis)**2/(2*(Lc**2)))

    """
    #Pf
    # background error Pf
    Pf = np.zeros((N,N))
    for k in range(en_member):
        e_f_k = (x_b_each[:,k]-x_b_mean[:]).reshape((N,1))
        Pf += e_f_k.dot(e_f_k.transpose())
    Pf /= (en_member-1)
    """    
    #PH
    PH = np.zeros((N,N))
    for k in range(en_member):
        c = (x_b_each[:,k]-x_b_mean[:]).reshape((N,1))
        d = (y_b[:,k]-y_b_mean[:]).reshape((1,N))
        PH += c.dot(d)
    PH /= (en_member-1)
    PH[:,:] = c_bloc[:,:]*PH[:,:]
    
    #HPH
    HPH = np.zeros((N,N))
    HPH = H.dot(PH)
    
    #Kalmen Gain K = PH*(HPH+R)^-1
    K_EnKF = PH.dot(np.linalg.inv(HPH+R))

    # analysis scheme
    x_a = x_b_each + K_EnKF.dot(y_o_en - y_b)
    
    #RTPP
    x_a = (1-alpha_en)*x_a + alpha_en*x_b_each
    x_a_en_save = np.vstack([x_a_en_save,x_a[np.newaxis,:,:]])    
    
    tt += 1  


x_a_mean = np.zeros((nT+1,N))
x_b_mean = np.zeros((nT+1,N))
for i in range(nT+1):
    for j in range(N):
        x_a_mean[i,j] = np.mean(x_a_en_save[i,j,:])
        x_b_mean[i,j] = np.mean(x_b_en_save[i,j,:])
        

# save EnKF background and analysis data(with inflation)

np.savetxt('x_b_EnKF.txt', x_b_mean)
np.savetxt('x_a_EnKF.txt', x_a_mean)
"""
np.savetxt('x_b_EnKF_obs1.txt', x_b_mean)
np.savetxt('x_a_EnKF_obs1.txt', x_a_mean)
"""
