"""
Plot the data assimilation results: pure localiz.
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
x_a_EnKF_save = np.genfromtxt('x_a_EnKF.txt')
x_a_EnKF04_save = np.genfromtxt('x_a_EnKF04.txt')
x_a_EnKF02_save = np.genfromtxt('x_a_EnKF02.txt')
x_a_EnKF01_save = np.genfromtxt('x_a_EnKF01.txt')
x_o_save= np.genfromtxt('x_o_save.txt')
x_no_da_save = np.genfromtxt('x_no_da_save.txt')

x = np.arange(nT+1) * dT
x_new = np.linspace(x.min(), x.max(), 1200)

# no DA RMSE
mean_noDA_err_rms = np.zeros(nT+1)
noDA_err_rms = np.sqrt((x_no_da_save-x_t_save)**2)
for i in range(nT+1):
    mean_noDA_err_rms[i] = np.mean(noDA_err_rms[i,:])
# EnKF RMSE
mean_EnKF_err_rms = np.zeros(nT+1)
EnKF_err_rms = np.sqrt((x_a_EnKF_save-x_t_save)**2)
for i in range(nT+1):
    mean_EnKF_err_rms[i] = np.mean(EnKF_err_rms[i,:])
# EnKF04 RMSE
mean_EnKF04_err_rms = np.zeros(nT+1)
EnKF04_err_rms = np.sqrt((x_a_EnKF04_save-x_t_save)**2)
for i in range(nT+1):
    mean_EnKF04_err_rms[i] = np.mean(EnKF04_err_rms[i,:])
# EnKF02 RMSE
mean_EnKF02_err_rms = np.zeros(nT+1)
EnKF02_err_rms = np.sqrt((x_a_EnKF02_save-x_t_save)**2)
for i in range(nT+1):
    mean_EnKF02_err_rms[i] = np.mean(EnKF02_err_rms[i,:])
# EnKF01 RMSE
mean_EnKF01_err_rms = np.zeros(nT+1)
EnKF01_err_rms = np.sqrt((x_a_EnKF01_save-x_t_save)**2)
for i in range(nT+1):
    mean_EnKF01_err_rms[i] = np.mean(EnKF01_err_rms[i,:])
    
#define spline
spl_noDA = make_interp_spline(x, mean_noDA_err_rms, k=7)
smooth_noda = spl_noDA(x_new)
spl_EnKF = make_interp_spline(x, mean_EnKF_err_rms, k=7)
smooth_EnKF = spl_EnKF(x_new)
spl_EnKF04 = make_interp_spline(x, mean_EnKF04_err_rms, k=7)
smooth_EnKF04 = spl_EnKF04(x_new)
spl_EnKF02 = make_interp_spline(x, mean_EnKF02_err_rms, k=7)
smooth_EnKF02 = spl_EnKF02(x_new)
spl_EnKF01 = make_interp_spline(x, mean_EnKF01_err_rms, k=7)
smooth_EnKF01 = spl_EnKF01(x_new)


fig, ax1= plt.subplots(1,1,sharex=True,sharey=False)
ax1.set_title(r'RMSE(Pure Localization; ensemble member:20)', size=20)
ax1.set_ylim(0,5.5)
ax1.set_xlim(0,10)
ax1.grid(color = 'gray', linestyle = '--', linewidth = 1)
ax1.plot(x_new, smooth_noda, 'g-' , label=r'$RMSE$_noDA')
ax1.plot(x_new, smooth_EnKF02, 'k-' , label=r'$RMSE$_EnKF(L=2)')
ax1.plot(x_new, smooth_EnKF, 'r-' , label=r'$RMSE$_EnKF(L=3)')
ax1.plot(x_new, smooth_EnKF04, 'b-' , label=r'$RMSE$_EnKF(L=4)')
ax1.plot(x_new, smooth_EnKF01, 'm-' , label=r'$RMSE$_EnKF(L=5)')
#ax1.plot(x_new, smooth_EnKF04, 'b-' , label=r'$RMSE$_EnKF(L=6)')
#ax1.plot(x_new, smooth_EnKF02, 'k-' , label=r'$RMSE$_EnKF(L=7)')
ax1.legend(loc='upper right', numpoints=1, prop={'size':18})
ax1.set_xlabel(r'$t$', size=18)
fig.subplots_adjust()
fig.show()
