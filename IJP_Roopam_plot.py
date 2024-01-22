## mechanical threshold model applied to SS316 SLM alloy with cellular sructure
import numpy as np
import matplotlib.pyplot as plt

########################
# Constants
C_Cr_in = 17
C_Mo_in = 2.5
C_Cr_wall = 18
C_Mo_in = 3
mu_0 = 71.49e9 #GPa
D_0 = 2.91e9 #GPa
T_0 = 204 #K
T = 300 #K
sigma_a = 40e6 #MPa
sigma_i = 450e6 #MPa
b = 2.49e-10 #Angstrom
g_0_i = 0.325
g_0_e = 1.6
eps_0_dot_i = 10e7 #s-1
eps_0_dot_e = 10e7 #s-1
eps_0_dot_s = 10e7 #s-1
p_i = 1/2
p_e = 2/3
q_i = 3/2
q_e = 1
h_0 = 4.1e9 #GPa
k= 1.38e-23
sigma_e_s_0 = 3.475e9 #GPa
g_0_e_s = 0.165
eps_dot = np.array([1e-3, 1e-2, 1e-1, 1, 10])

#######################
mu = mu_0 -D_0/(np.exp(T_0/T)-1)
s_i = np.power(1-np.power((k*T/(mu*b*b*b*g_0_i))*np.log(eps_0_dot_i/eps_dot),1/q_i),1/p_i)
s_e = np.power(1-np.power((k*T/(mu*b*b*b*g_0_e))*np.log(eps_0_dot_e/eps_dot),1/q_e),1/p_e)
# sigma_f = sigma_a +(mu/mu_0)(s_i*sigma_i+s_e*sigma_e)

figure1, ax = plt.subplots(1,1,figsize=(8,6))#, layout='constrained')
ax.plot(eps_dot,s_i,'+b')
ax.plot(eps_dot,s_e,'*r')
ax.set_xscale('log')
ax.set_xlabel('$\\dot{\\epsilon}$ (s$^{-1}$)',fontsize=12)
ax.set_ylabel('$s_{i}, s_{e}$',fontsize=14)
ax.set(xlim = [1e-4, 10],ylim = [0, 1])    
ax.legend(['$s_{i}$', '$s_{e}$'],loc='lower right',fontsize=12)
ax.tick_params(axis='both', which='both', direction='in',labelsize=14,top=True, right=True)
ax.set_xticks([1e-3, 1e-2, 1e-1, 1,10])
ax.set_xticklabels(['10$^{-3}$', '10$^{-2}$', '10$^{-1}$', '1', '10'])
ax.set_yticks([0.25, 0.5, 0.75, 1])
ax.set_yticklabels(['0.25', '0.50', '0.75', '1'])
figure1.savefig(r'C:\Users\MKY\OneDrive\PythonWork\utilities\figureReviewer.tif',dpi = 400,)