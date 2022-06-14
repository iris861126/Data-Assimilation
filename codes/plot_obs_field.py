"""
Plot the data assimilation results: obs. field
    
Read:
  x_t.txt
  x_b.txt
  x_a.txt
  x_o_save.txt
"""
import numpy as np
from settings import *
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.interpolate import make_interp_spline, BSpline

# load data
x_t_save = np.genfromtxt('x_t.txt')
x_a_oi_save = np.genfromtxt('x_a_oi.txt')
x_a_oi_obs_save = np.genfromtxt('x_a_oi_obs1.txt')
x_a_EKF_infla_save = np.genfromtxt('x_a_EKF_infla.txt')
x_a_EKF_obs_infla_save = np.genfromtxt('x_a_EKF_infla_obs1.txt')
x_a_in3d_save = np.genfromtxt('x_a_in3d.txt')
x_a_in3d_obs_save = np.genfromtxt('x_a_in3d_obs1.txt')
x_a_EnKF_save = np.genfromtxt('x_a_EnKF.txt')
x_a_EnKF_obs_save = np.genfromtxt('x_a_EnKF_obs1.txt')
x_o_save= np.genfromtxt('x_o_save.txt')
x_no_da_save = np.genfromtxt('x_no_da_save.txt')

Pb = np.genfromtxt('Pb.txt')

x = np.arange(nT+1) * dT
x_new = np.linspace(x.min(), x.max(), 1200)

# no DA RMSE
mean_noDA_err_rms = np.zeros(nT+1)
noDA_err_rms = np.sqrt((x_no_da_save-x_t_save)**2)
for i in range(nT+1):
    mean_noDA_err_rms[i] = np.mean(noDA_err_rms[i,:])

# OI RMSE
mean_OI_err_rms = np.zeros(nT+1)
OI_err_rms = np.sqrt((x_a_oi_save-x_t_save)**2)
for i in range(nT+1):
    mean_OI_err_rms[i] = np.mean(OI_err_rms[i,:])   
# OI_obs RMSE
mean_OI_obs_err_rms = np.zeros(nT+1)
OI_obs_err_rms = np.sqrt((x_a_oi_obs_save-x_t_save)**2)
for i in range(nT+1):
    mean_OI_obs_err_rms[i] = np.mean(OI_obs_err_rms[i,:])
# incremental 3dVar RMSE
mean_in3d_err_rms = np.zeros(nT+1)
in3d_err_rms = np.sqrt((x_a_in3d_save-x_t_save)**2)
for i in range(nT+1):
    mean_in3d_err_rms[i] = np.mean(in3d_err_rms[i,:])
# incremental 3dVar_obs RMSE
mean_in3d_obs_err_rms = np.zeros(nT+1)
in3d_obs_err_rms = np.sqrt((x_a_in3d_obs_save-x_t_save)**2)
for i in range(nT+1):
    mean_in3d_obs_err_rms[i] = np.mean(in3d_obs_err_rms[i,:])

# EKF RMSE    
mean_EKF_infla_err_rms = np.zeros(nT+1)
EKF_infla_err_rms = np.sqrt((x_a_EKF_infla_save-x_t_save)**2)
for i in range(nT+1):
    mean_EKF_infla_err_rms[i] = np.mean(EKF_infla_err_rms[i,:])    
# EKF_obs RMSE
mean_EKF_obs_infla_err_rms = np.zeros(nT+1)
EKF_obs_infla_err_rms = np.sqrt((x_a_EKF_obs_infla_save-x_t_save)**2)
for i in range(nT+1):
    mean_EKF_obs_infla_err_rms[i] = np.mean(EKF_obs_infla_err_rms[i,:])
# EnKF RMSE
mean_EnKF_err_rms = np.zeros(nT+1)
EnKF_err_rms = np.sqrt((x_a_EnKF_save-x_t_save)**2)
for i in range(nT+1):
    mean_EnKF_err_rms[i] = np.mean(EnKF_err_rms[i,:])
# EnKF_obs RMSE
mean_EnKF_obs_err_rms = np.zeros(nT+1)
EnKF_obs_err_rms = np.sqrt((x_a_EnKF_obs_save-x_t_save)**2)
for i in range(nT+1):
    mean_EnKF_obs_err_rms[i] = np.mean(EnKF_obs_err_rms[i,:])



#define spline
spl_noDA = make_interp_spline(x, mean_noDA_err_rms, k=7)
smooth_noda = spl_noDA(x_new)
#-----------------------
spl_OI = make_interp_spline(x, mean_OI_err_rms, k=7)
smooth_OI = spl_OI(x_new)

spl_OI_obs = make_interp_spline(x, mean_OI_obs_err_rms, k=7)
smooth_OI_obs = spl_OI_obs(x_new)

spl_in3d = make_interp_spline(x, mean_in3d_err_rms, k=7)
smooth_in3d = spl_in3d(x_new)

spl_in3d_obs = make_interp_spline(x, mean_in3d_obs_err_rms, k=7)
smooth_in3d_obs = spl_in3d_obs(x_new)
#-----------------------
spl_EKF_infla = make_interp_spline(x, mean_EKF_infla_err_rms, k=7)
smooth_EKF_infla = spl_EKF_infla(x_new)

spl_EKF_obs_infla = make_interp_spline(x, mean_EKF_obs_infla_err_rms, k=7)
smooth_EKF_obs_infla = spl_EKF_obs_infla(x_new)

spl_EnKF = make_interp_spline(x, mean_EnKF_err_rms, k=7)
smooth_EnKF = spl_EnKF(x_new)

spl_EnKF_obs = make_interp_spline(x, mean_EnKF_obs_err_rms, k=7)
smooth_EnKF_obs = spl_EnKF_obs(x_new)
#-----------------------
# Plot time series of grid point avg.
fig, ax1= plt.subplots(1,1,sharex=True,sharey=False)
#ax1.set_title(r'RMSE(test:[a]; OI & 3DVAR)', size=20)
ax1.set_ylim(0,5.5)
ax1.set_xlim(0,10)
ax1.grid(color = 'gray', linestyle = '--', linewidth = 1)
ax1.plot(x_new, smooth_noda, 'g-' , label=r'$RMSE$_noDA')
"""
#----------
#OI, 3DVAR
ax1.set_title(r'RMSE(test:[b]; OI & 3DVAR)', size=20)
ax1.plot(x_new, smooth_in3d,color='lightcoral', linewidth=3,linestyle='-', label=r'$RMSE$_3DVAR')
ax1.plot(x_new, smooth_in3d_obs,color='#D2A2CC', linewidth=3,linestyle='-', label=r'$RMSE$_3DVAR(modi. obs.)')
ax1.plot(x_new, smooth_OI,color='#4b0082', linestyle='--', label=r'$RMSE$_OI')
ax1.plot(x_new, smooth_OI_obs,color='#3D7878', linestyle='--', label=r'$RMSE$_OI(modi. obs.)')
#----------

"""
#----------
#EnKF, EKF
ax1.set_title(r'RMSE(test:[b]; EKF & EnKF)', size=20)
ax1.plot(x_new, smooth_EnKF_obs, 'm-' , label=r'$RMSE$_EnKF(modi. obs.)')
ax1.plot(x_new, smooth_EnKF, color='#613030', linestyle='-', label=r'$RMSE$_EnKF')
ax1.plot(x_new, smooth_EKF_obs_infla, color='#CF9E9E',linestyle='-', label=r'$RMSE$_EKF(modi. obs.)')
ax1.plot(x_new, smooth_EKF_infla, color='#9F4D95',linestyle='-', label=r'$RMSE$_EKF')
#----------

ax1.legend(loc='upper right', numpoints=1, prop={'size':18})
ax1.set_xlabel(r'$t$', size=18)
fig.subplots_adjust()
fig.show()
