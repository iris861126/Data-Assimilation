"""
The data assimilation system:OI
Load:
  x_a_init.txt
Save:
  x_b_oi.txt
  x_a_oi.txt
"""
import numpy as np
from scipy.integrate import ode
import lorenz96
from settings import *
import math
import seaborn as sns

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


# load truth
x_t_save = np.genfromtxt('x_t.txt')

# load H
#H = np.genfromtxt('H2.txt')
H = np.eye(N)

# calculate for obs. error covariance(R)

R = np.zeros((N,N))
for i in range(N):
#    R[i,i] = np.var(y_o_save[:,i])
     R[i,i] = (std_obs_err)**2.
Pb = np.zeros((N,N))
Pb[:,:] = R[:,:]    #let Pb = R 
"""
# density obs. err 
R = np.zeros((19,19))
for i in range(19):
     R[i,i] = (std_obs_err)**2.

Pb = np.zeros((N,N))
Pb = np.transpose(H).dot(R.dot(H))    #let Pb = R 
"""
# NMC process
iter = 0
while True:
    tt = 1    
    # initial x_b: no values at the initial time (assign NaN)
    x_b_save = np.full((1,N), np.nan, dtype='f8')
    # initial x_a: from x_a_init
    x_a_save = np.array([x_a_init])
    while tt <= nT:
        print('DA process for NMC...')
        tts = tt - 1
        Ts = tts * dT  # forecast start time
        Ta = tt  * dT  # forecast end time (DA analysis time)
#        print('Cycle =', tt, ', Ts =', round(Ts, 10), ', Ta =', round(Ta, 10))
        
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
        
        # obs. err
        #obs_err = y_o-x_t
        
        # innovation
        y_b = np.dot(H, x_b)
        d = y_o - y_b
        
        #Kalmen Gain
        #K_OI = np.dot((np.dot(Pb,H.transpose())),(np.linalg.inv((np.dot(H,np.dot(Pb,H.transpose())))+R)))
        K_OI = Pb.dot(H.transpose()).dot(np.linalg.inv(H.dot(Pb.dot(H.transpose()))+R))
    
        # analysis scheme
        x_a = x_b + np.dot(K_OI,d)

        x_a_save = np.vstack([x_a_save, x_a.transpose()])
        tt += 1

    if(iter < 5):
        iter += 1
    else: 
        print("Pb= ",Pb)
        break    

    #the NMC method
    print('NMC method...')
    xf_t1 = np.zeros((N))
    xf_t2 = np.zeros((N))
    tt_NMC = 1    
    while tt_NMC <= nT_NMC:
        tts = tt_NMC - 1
        tts_NMC = factor_NMC*tts
        ttf = tt_NMC + 1
        Ts_NMC = tts * dT_NMC  # forecast start time
        Tf_NMC = ttf * dT_NMC  # forecast end time
#        print('Cycle =', tt_NMC, ', tts =', round(tts, 10), ', Ts_NMC =', round(Ts_NMC, 10) , ', Tf_NMC =', round(Tf_NMC, 10))
    
        #-------------------
        # extended forecast 
        #-------------------
    
        solver = ode(lorenz96.f).set_integrator('dopri5')
        solver.set_initial_value(x_a_save[tts_NMC], Ts_NMC).set_f_params(F)
        solver.integrate(Tf_NMC-dT_NMC)
        xf_t1 = np.vstack([xf_t1, [solver.y]]) 
        solver.integrate(Tf_NMC)
        xf_t2 = np.vstack([xf_t2, [solver.y]])
    
        tt_NMC += 1
        
    xf_t1 = np.delete(xf_t1, (0, 1), 0) #saving row2 to row50
    xf_t2 = np.delete(xf_t2, (0,num_NMC), 0) #saving row1 to row49
    t2_minus_t1 = xf_t1 - xf_t2
    t2_minus_t1_T = np.transpose(t2_minus_t1)     
    Pb = alpha*(np.dot(t2_minus_t1_Ｔ,t2_minus_t1))/(num_NMC-1)   

    # average for BEC
    Pbsum = np.zeros(N)
    times = np.zeros(N)
    for m in range(N):
        x = round((m/math.sqrt(2)),3)
        for i in range(N):
            for j in range(N):
                if(round(getDis(i,j),3) == x):   #using distance to classify them
                    Pbsum[m] = Pbsum[m]+Pb[i,j]
                    times[m] += 1
    #for row1
    for m in range(N):
        x = round((m/math.sqrt(2)),3)
        for i in range(21):
            j=0
            if(round(getDis(i,j),3) == x):
                Pb[i,j] = Pbsum[m]/times[m]
    for i in range(1,21):
        pt = N+(-i)
        Pb[pt,j]=Pb[i,j]
    Pb = Pb.transpose()   
    
    #for row2 t0 row40
    for i in range(1,N):
        for j in range(N):
            if((j-1)<0):
                j2 = N+(j-1)
                Pb[i,j] = Pb[i-1,j2]
            else:
                Pb[i,j] = Pb[i-1,j-1]
np.savetxt('Pb.txt', Pb)

# DA process:OI
tt = 1    
# initial x_b: no values at the initial time (assign NaN)
x_b_save = np.full((1,N), np.nan, dtype='f8')
# initial x_a: from x_a_init
x_a_save = np.array([x_a_init])
while tt <= nT:
    print('DA process...')
    tts = tt - 1
    Ts = tts * dT  # forecast start time
    Ta = tt  * dT  # forecast end time (DA analysis time)
#        print('Cycle =', tt, ', Ts =', round(Ts, 10), ', Ta =', round(Ta, 10))
    
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
    
    # obs. err
    #obs_err = y_o-x_t
    
    # innovation
    y_b = np.dot(H, x_b)
    d = y_o - y_b
    
    #Kalmen Gain
    K_OI = np.dot((np.dot(Pb,np.transpose(H))),(np.linalg.inv((np.dot(H,np.dot(Pb,np.transpose(H))))+R)))

    # analysis scheme
    x_a = x_b + np.dot(K_OI,d)

    x_a_save = np.vstack([x_a_save, x_a.transpose()])
    tt += 1

# save OI background and analysis data
np.savetxt('x_b_oi.txt', x_b_save)
np.savetxt('x_a_oi.txt', x_a_save)
"""
# save OI background and analysis data
np.savetxt('x_b_oi_obs1.txt', x_b_save)
np.savetxt('x_a_oi_obs1.txt', x_a_save)
"""
#%%
sns.heatmap(Pb, cmap='RdYlGn', xticklabels=5, yticklabels=5)
