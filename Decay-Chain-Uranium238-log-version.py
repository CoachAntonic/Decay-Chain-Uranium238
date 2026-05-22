import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def chain_decay(y, time_s) :
    #Decay constant of element in the chain decay of the Uranium 238 in s^-1
    lambda_U238 = np.log(2) / 1.410e17 
    lambda_Th234 = np.log(2) / 2082240 
    lambda_Pa234 = np.log(2) / 70.2    
    lambda_U234 = np.log(2) / 7.748e12 
    lambda_Th230 = np.log(2) / 2.349e12 
    lambda_Ra226 = np.log(2) / 5.049e10 
    lambda_Rn222 = np.log(2) / 330048 
    lambda_Po218 = np.log(2) / 183  
    lambda_Pb214 = np.log(2) / 1608 
    lambda_Bi214 = np.log(2) / 1182 
    lambda_Po214 = np.log(2) / 0.000164
    lambda_Pb210 = np.log(2) / 7.034e8 
    lambda_Bi210 = np.log(2) / 433094 
    lambda_Po210 = np.log(2) / 11959000 
    #Pb 216 is stable 
    (U_238, Th_234, Pa_234, 
     U_234, Th_230, Ra_226, 
     Rn_222, Po_218, Pb_214, 
     Bi_214, Po_214, Pb_210,
     Bi_210, Po_210, Pb_206) = y
    
    dydt = [-lambda_U238 * U_238, lambda_U238 * U_238 - lambda_Th234 * Th_234,
     lambda_Th234 * Th_234 - lambda_Pa234 * Pa_234, lambda_Pa234 * Pa_234 - lambda_U234 * U_234,
     lambda_U234 * U_234 - lambda_Th230 * Th_230, lambda_Th230 * Th_230 - lambda_Ra226 * Ra_226,
     lambda_Ra226 * Ra_226 - lambda_Rn222 * Rn_222, lambda_Rn222 * Rn_222 - lambda_Po218 * Po_218,
     lambda_Po218 * Po_218 - lambda_Pb214 * Pb_214, lambda_Pb214 * Pb_214 - lambda_Bi214 * Bi_214,
     lambda_Bi214 * Bi_214 - lambda_Po214 * Po_214, lambda_Po214 * Po_214 - lambda_Pb210 * Pb_210,
     lambda_Pb210 * Pb_210 - lambda_Bi210 * Bi_210, lambda_Bi210 * Bi_210 - lambda_Po210 * Po_210, 
     lambda_Po210 * Po_210
]
    return dydt
 
#The value of these variables are changeable
half_life_U_238 = 1.410e17
time_s = np.linspace(0, 20*half_life_U_238, 1000)
Atoms_U_238_at_time_0 = 10000

#For this exemple we assume that the number of atoms for all the isotopes exepts the urabium 238 are 0 at time 0
y0 = [Atoms_U_238_at_time_0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

#Now we calculate our solution
solution = odeint(chain_decay, y0, time_s)

#We do 8 graphs for 8 groups of elements with an half life in the same range
plt.style.use('dark_background')
fig, axes = plt.subplots(4, 2, figsize=(22, 15))

#Graph for the uranium 238 and plomb 206
axes[0, 0].plot(time_s, solution[:,0], 'r-', label='Uranium 238')
axes[0, 0].plot(time_s, solution[:,14], 'w-', label='Lead 206')
axes[0,0].set_xscale('log')
axes[0,0].set_yscale('log')
axes[0, 0].set_xlabel('Time (s)')
axes[0, 0].set_ylabel('N (atoms)')
axes[0,0].legend()

#Graph for the Uranium 234, Thorium 230
axes[0, 1].plot(time_s, solution[:,3], 'g-', label='Uranium 234')
axes[0, 1].plot(time_s, solution[:,4], 'w-', label='Thorium 230')
axes[0, 1].set_xscale('log')
axes[0, 1].set_yscale('log')
axes[0, 1].set_xlabel('Time (s)')
axes[0, 1].set_ylabel('N (atoms)')
axes[0,1].legend()

#Graph for the radium 226 
axes[1, 0].plot(time_s, solution[:,5], 'c-', label='Radium 226')
axes[1, 0].set_xscale('log')
axes[1, 0].set_yscale('log')
axes[1, 0].set_xlabel('Time (s)')
axes[1, 0].set_ylabel('N (atoms)')
axes[1,0].legend()

#Graph for the Lead 210
axes[1, 1].plot(time_s, solution[:,11], 'r-', label='Lead 210')
axes[1, 1].set_xscale('log')
axes[1, 1].set_yscale('log')
axes[1, 1].set_xlabel('Time (s)')
axes[1, 1].set_ylabel('N (atoms)')
axes[1,1].legend()

#Graph for the Polonium 210 and Thorium 234
axes[2, 0].plot(time_s, solution[:,13], 'b-', label='Polonium 210')
axes[2, 0].plot(time_s, solution[:,1], 'w-', label='Thorium 234')
axes[2, 0].set_xscale('log')
axes[2, 0].set_yscale('log')
axes[2, 0].set_xlabel('Time (s)')
axes[2, 0].set_ylabel('N (atoms)')
axes[2,0].legend()

#Graph for the Bismuth 210 and Radon 222
axes[2, 1].plot(time_s, solution[:,12], 'r-', label='Bismuth 210')
axes[2, 1].plot(time_s, solution[:,6], 'w-', label='Radon 222')
axes[2, 1].set_xscale('log')
axes[2, 1].set_yscale('log')
axes[2, 1].set_xlabel('Time (s)')
axes[2, 1].set_ylabel('N (atoms)')
axes[2,1].legend()

#Graph for the Lead 214, Polonium 218, Bismuth 214 and Protactinium 234
axes[3, 0].plot(time_s, solution[:,8], 'r-', label='Lead 214')
axes[3, 0].plot(time_s, solution[:,7], 'w-', label='Polonium 218')
axes[3, 0].plot(time_s, solution[:,9], 'g-', label='Bismuth 214')
axes[3, 0].plot(time_s, solution[:,2], 'c-', label='Protactinium 234')
axes[3, 0].set_xscale('log')
axes[3, 0].set_yscale('log')
axes[3, 0].set_xlabel('Time (s)')
axes[3, 0].set_ylabel('N (atoms)')
axes[3,0].legend()

#Graph for the Polonium 214
axes[3, 1].plot(time_s, solution[:,10], 'r-', label='Polonium 214')
axes[3, 1].set_xscale('log')
axes[3, 1].set_yscale('log')
axes[3, 1].set_xlabel('Time (s)')
axes[3, 1].set_ylabel('N (atoms)')
axes[3,1].legend()

plt.suptitle('Complete decay chain of Uranium 238', fontsize=40)
plt.tight_layout()
plt.show()


    
    
    


