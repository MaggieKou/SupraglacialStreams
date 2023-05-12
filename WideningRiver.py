# -*- coding: utf-8 -*-
"""
Created on Tue May  9 10:54:56 2023

@author: maggi
"""
import numpy as np
from matplotlib import pyplot as plt 

#numerical method to implement on equation 21 from Moulins.pdf by Neufeld and Devauchelle
#given that the widening river is an ODE of S: dS/dt = -5/2 S^(12/5)

n = 10000 #num of steps 
#assuming t0 = 0
tf = 100
del_t = tf/n #step size, assuming uniform
S_o = 5 #initial slope 

t = np.linspace(0, tf, n)
S = np.zeros([n])
S[0] = S_o

def fwd_euler(t, S):
    for i in range (1, n):
        S[i] = S[i-1] - del_t*(5/2*np.sign(S[i-1])*np.abs(S[i-1])**(12/5))
    return S

#plot ananlytic solution 
S_analytic = np.zeros([n])
S_an2 = np.zeros([n])
for i in range (0,n):
    S_analytic[i] = (S_o**(7/5)-7/2*t[i])**(5/7)
    S_an2[i] = (7/2*t[i] + S_o**(-7/5))**(-5/7)

S_approx = fwd_euler(t, S)
print (t , "\n" , S)

fig, ax = plt.subplots()
plt_analytic = ax.plot(t, S_analytic, '-', label="Analytic Solution")
plt_an2 =  ax.plot(t, S_an2, '-', label="Analytic Solution2")
plt_approx = ax.plot(t, S_approx, '--', label = "Numerical approximation with Forward Euler")
ax.legend( ["Analytic Solution (22) in paper","Alternate Analytic Solution","Numerical approximation with Forward Euler Method"])

plt.xlabel("time [dimensionless]")
plt.ylabel("slope [dimensionless]")
plt.title("Euler Fwd on Widening River Case")
plt.show()


    
     
