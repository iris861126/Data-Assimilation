"""
General settings
"""
import numpy as np

#===================
# global parameters
#===================
N = 40                      # number of grid point
F = 8.                      # forcing term
Tmax = 10.                  # time length of the experiment (50 days)
dT = 0.05                   # DA cycle length (1=5days -> 6hrs=0.05)
nT = int(Tmax/dT)           # number of cycles
#setting for obs.
mean_obs_err = 0.
std_obs_err = 2.

#=================
# setting for NMC
#=================
factor_NMC = 4
dT_NMC = factor_NMC*dT      #0.2(24hr)
nT_NMC = int(Tmax/dT_NMC) 
alpha = 0.4                 #0.3, 0.65
num_NMC = int(nT/factor_NMC)


#=================
# setting for EKF
#=================
num_TLM= 6                  #seperate 1 time step into [num_TLM] intervals (eg. 6)
dT_TLM = dT/num_TLM         #for num_TLM=5, dT_TLM=0.01(1.2hr)
nT_TLM = int(dT/dT_TLM) 
ro_EKF = 1.5                #for multiplicative infla.
add_factor_infla = 0.       #for additive infla.

#===================
# setting for 4dVar
#===================
factor_4dVar = 4
dT_4d = dT/factor_4dVar
nT_4d = int(Tmax/dT_4d)
#for 4dTLM
num_4dTLM= 8                  
dT_4dTLM = 2*dT/num_4dTLM     
nT_4dTLM = 2*int(dT/dT_4dTLM)
num_L_4dTLM = 5
dT_L_4dTLM = dT_4dTLM/num_L_4dTLM 
nT_L_4dTLM = int(dT_4dTLM/dT_L_4dTLM)

#==================
# setting for EnKF
#==================
en_member = 100 
mean_en_initial_err = 0.         # ensemble mean = the initial condi. of other DA method
std_en_initial_err = [2.]*N
alpha_en = 0.5                   #0.2,0.5,0.7
Lc = 7






