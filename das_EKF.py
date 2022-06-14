"""
The data assimilation system:EKF
Load:
  x_a_init.txt
Save:
  x_b_EKF.txt
  x_a_EKF.txt
"""
import numpy as np
from scipy.integrate import ode
import lorenz96
from settings import *
import math

def getDis(pointX,pointY):
    a=39. - 1.          #lineY2-lineY1
    b=1. - 39.          #lineX1-lineX2
    c=39.*1. - 1.*39. #lineX2*lineY1-lineX1*lineY2
    dis=(math.fabs(a*pointX+b*pointY+c))/(math.pow(a*a+b*b,0.5))
    return dis

# load initial condition
x_a_init = np.genfromtxt('x_a_init.txt')

# load observations
y_o_save = np.genfromtxt('x_o_save.txt')

# initial x_b: no values at the initial time (assign NaN)
x_b_save = np.full((1,N), np.nan, dtype='f8')

# initial x_a: from x_a_init
x_a_save = np.array([x_a_init])

# load truth
x_t_save = np.genfromtxt('x_t.txt')

#load initial BEC
Pb_EKF = np.genfromtxt('Pb.txt')    #initial BEC

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
#DA process: EKF
tt = 1    
Pa_EKF = Pb_EKF
while tt <= nT:
    print('EKF DA process...')
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
    
    #TLM
    print('TLM process...')
    x_a_state = x_a_save[tts]
    Ts_TLM=Ts
    Ta_TLM = Ts_TLM + dT_TLM  # TLM forecast end time 
    print(' Ts_TLM =', round(Ts_TLM, 10), ', Ta_TLM =', round(Ta_TLM,10))
    tt_TLM = 1
    L_final = np.zeros((N,N,num_TLM))
    L_old = np.zeros((N,N))
    while tt_TLM <= nT_TLM:      
        tts_TLM = tt_TLM-1
        L = np.zeros((N,N))
        dt_F = np.zeros((N,N))
        for i in range(N): 
            if(i-1<0):
                a1 = N+(i-1)
            else:
                a1 = i-1
            if(i-2<0):
                a2 = N+(i-2)
            else:
                a2 = i-2
            if(i+1>=N):
                b1 = -N+(i+1)
            else:
                b1 = i+1
            #print("a1 a2 b1=",a1,a2,b1)
            
            dt_F[i] = [-dT_TLM, dT_TLM*x_a_state[a1]]+[0.]*36+[-dT_TLM*x_a_state[a1], dT_TLM*(x_a_state[b1]-x_a_state[a2])]
            dt_F[i] = np.roll(dt_F[i],i) 
            #print('i= ',i)
        L = np.eye(N) + dt_F
        L_old = L
        Ts_TLM = Ta_TLM
        Ta_TLM = Ts_TLM + dT_TLM
        if(tt_TLM < num_TLM):
            solver = ode(lorenz96.f).set_integrator('dopri5')
            solver.set_initial_value(x_a_state, Ts_TLM).set_f_params(F)
            print(' Ts_TLM =', round(Ts_TLM,10), ', Ta_TLM =', round(Ta_TLM,10))
            x_a_state = solver.integrate(Ta_TLM)
            L_final[:,:,tts_TLM] = L
            tt_TLM += 1
        else:
            L_final[:,:,tts_TLM] = L
            break

    M = L_final[:,:,(num_TLM-1)]  #start from the last L
    for j in np.arange(num_TLM-2,-1,-1):
        #M = np.dot(M,L_final[:,:,j])
        M = L_final[:,:,j].dot(M)

    Pb_EKF = np.dot(np.dot(M,Pa_EKF),np.transpose(M))                 

    #--------------
    # analysis step
    #--------------
        
    #truth
    x_t = x_t_save[tt].transpose()
    # background
    x_b = x_b_save[tt].transpose()
    
    # observation
    y_o = y_o_save[tt].transpose()
    
    # obs. err
    #obs_err = y_o-x_t
    
    # innovation
    y_b = np.dot(H, x_b)
    d = y_o - y_b
    
    #Kalmen Gain
    K_EKF = np.dot((np.dot(Pb_EKF,np.transpose(H))),(np.linalg.inv((np.dot(H,np.dot(Pb_EKF,np.transpose(H))))+R)))

    # analysis scheme
    x_a = x_b + np.dot(K_EKF,d)
    x_a_save = np.vstack([x_a_save, x_a.transpose()])
    
    Pa_EKF = np.dot((np.eye(N) - np.dot(K_EKF,H)),Pb_EKF)
    #BEC inflation
    Pa_EKF = ro_EKF*Pa_EKF #multiplicative infla.
    tt += 1  


if(ro_EKF == 1.): 
    # save EKF background and analysis data(no inflation)
    np.savetxt('x_b_EKF.txt', x_b_save)
    np.savetxt('x_a_EKF.txt', x_a_save)
    print("No infla. data is saved.")
elif(ro_EKF != 1.):
    # save EKF background and analysis data(with inflation)
    np.savetxt('x_b_EKF_infla.txt', x_b_save)
    np.savetxt('x_a_EKF_infla.txt', x_a_save)
    print("Infla. data is saved.")
"""
if(ro_EKF == 1.): 
    # save EKF background and analysis data(no inflation)
    np.savetxt('x_b_EKF_obs1.txt', x_b_save)
    np.savetxt('x_a_EKF_obs1.txt', x_a_save)
    print("No infla. data is saved.")
elif(ro_EKF != 1.):
    # save EKF background and analysis data(with inflation)
    np.savetxt('x_b_EKF_infla_obs1.txt', x_b_save)
    np.savetxt('x_a_EKF_infla_obs1.txt', x_a_save)
    print("Infla. data is saved.")
"""