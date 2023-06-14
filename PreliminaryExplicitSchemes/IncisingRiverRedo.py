# -*- coding: utf-8 -*-
"""
Created on Wed May 31 12:08:27 2023

@author: maggi
"""

import numpy as np
from matplotlib import pyplot as plt
"Set the solver parameters here"

#TIMESTEPS
to = 0 #starting at the intial condition
tf = 10
n_t = 1000 #number of timesteps 
del_t = (tf-to)/n_t
#ho = np.linspace(to, tf, n_t) #initiate initial conditions 
t = np.linspace(to, tf, n_t)

#SPATIAL STEPS
xo = 0 #starting at the intial condition
xf = 0.5
n_x = 10 #number of spatial steps
del_x = (xf-xo)/n_x
#ho = np.linspace(to, tf, n_t) #initiate initial conditions 
x = np.linspace(xo, xf, n_x)

#setting the parameter space
h = np.zeros([n_x, n_t])
S = np.zeros([n_x, n_t])

#solution convergence
ep_in = 0.0001 #parameter for initial boundary condition solution
ep_S = 0.0001 #parameter for the convergence of the final S 


#set up the solver for the boundary conditions on both h and S
R = 0.5 #surface stream


#INTIAL CONDITIONS 
#polynomial example to the 3rd degree
a = 0
b = -3
c = 0.1
d = 100
for i in range (0, n_x): #is this not being initialized ? 
    h[i][0] = a*x[i]**3 + b*x[i]**2 + c*x[i] + d # like a linear down sloping from the initial position
    S[i][0] = (h[i][0] - h[i-1][0])/del_x #slope at (0,0) is not defined through this 
S[0][0] = (h[1][0] - h[0][0])/del_x #Initialize first point

#BOUNDAY CONDITIONS

"constants that set initial conditions"
F = 1373 #W/m^2
rho = 0.001 #density of freshwaterkg/m^3
L = 3*10**5 #latent heat of the ice
Q = F/(rho*L) #slope 

So = S[0][0]
St = np.zeros([n_t])
#sets initial time derivative of the slope 
Sto = -(5/2*(Q**(1/5)*np.sign(So)*np.power(abs(So), 7/5)*np.sqrt(1 + So**2))/(R))
St[0] = Sto #this can be true for all x since the flux is constant across all the points along the river

 
for i in range (0, n_x):   
    for j in range (1, n_t):
        Snw = St[j-1] #new slope time derivative, to compare against
        St[j] = St[j-1] + 7/5*del_t*(np.power(St[j-1],2))/S[i][j-1] #solves for the derivative at this stage from previous "S" value
        S[i][j] = S[i][j-1] + del_t*St[j]#new slope at this point from previous time step 
        ht = -np.power(abs(S[i][j]), 6/5)/(np.sqrt(1 + np.power(abs(S[i][j]), 2))) #calculates time derivative
        h[i][j] = h[i][j-1] + del_t*ht #uses PDE as a corrector

for i in range (1, n_x):   
    for j in range (1, n_t):
        S_temp = (h[i][j] - h[i-1][j])/(del_x) #slope of the height, should always be negative
        
        while (S_temp - S[i][j] > ep_in):
            S[i][j] = S_temp #slope set to new iteration
            ht = -np.power(abs(S_temp), 6/5)/(np.sqrt(1 + np.power(abs(S_temp), 2))) #calculates time derivative
            h[i][j] = h[i][j-1] + del_t*ht #uses PDE as a corrector
            S_temp = (h[i][j] - h[i-1][j])/(del_x) #slope of the height, should always be negative
            print (S[i][j], S_temp)
        S[i][j] = S_temp

for i in range (1, n_x):
    for j in range (0, n_t):
        h[i][j] = h[i-1][j] + del_x*S[i-1][j]
        
    #this scheme is working, try receeding case

h_ints = np.zeros([5,len(x)]) 
time_label = 5*["a"]
for i in range (0, n_x):
    h_ints[0][i] = h[i][0]
for i in range (0, n_x): 
    for j in range (1, 5):
        h_ints[j][i] = h[i][int(n_t/(5 - j))-1] #samples the discrete times
        time_label[j] = "t = " + str(round(t[int(n_t/(5 - j))-1],2))
    #print(h_ints)
fig,ax = plt.subplots()

for plt_time in range (0, 5):
    h_plt = ax.plot(x[:], h_ints[plt_time][:]) #all the times plotted
  
ax.legend(time_label) #there may be an issue with global permissions for this to show up 

plt.xlabel("x [dimensionless]")
plt.ylabel("h [dimensionless]")
plt.title("Solving Scheme on h for different times")
plt.show()
   

