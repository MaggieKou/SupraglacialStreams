#Const Flux Assumption (i.e. in both space and time)

import numpy as np
import math
from matplotlib import pyplot as plt 

#ssuming constant slope case
F = 1373 #W/m^2
rho = 0.001 #density of freshwaterkg/m^3
L = 3*10**5 #latent heat of the ice
h_t = F/(rho*L) #slope 


#with scheme, solve after linearization
t_f = 2
t_0 = 0
n_t = 1000
del_t = (t_f - t_0)/n_t

#spatial steps
x_0 = 0
x_f = 5
n_x = 100
del_x = (x_f - x_0)/n_x

#initial conditions - assuming a linear profile of the form h = h_top - a*x
h_top = 1000
a = 0.1 #slope

h_in = h_top -a*np.linspace(x_0,x_f, n_x) #linear initial distribution
#solutions 
t = np.linspace(t_0,t_f, n_t) 
x = np.linspace(x_0,x_f, n_x) 

h_xt =np.zeros([n_t,n_x])
h_xt[0] = h_in #initial distribution

h_x = np.zeros([4,n_t]) #h_x at x = x_o for all t, need another relation to forward step this 

h_x_t = np.zeros([n_t-1]) #derivative of h_x wrt to t

xi = np.zeros([4,n_t]) 

xi[0][0] = -np.power(0.01,-1/5)
xi[1][0] = -np.power(0.1,-1/5)
xi[2][0] = -np.power(1,-1/5)
xi[3][0] = -np.power(2,-1/5)


h_x[0][0] = -0.01
h_x[1][0] = -0.1
h_x[2][0] = -1
h_x[3][0] = -2

for i in range (1, n_t):
    for k in range (0,4):
        xi[k][i]= xi[k][i-1] + del_t*(h_t/5 -1/2*(xi[k][i-1])**(-6)) #gives xi which is h_x^(-1/5)
        h_x[k][i]= np.power(xi[k][i],-5) #construct how slope varies with time, slope is spatially invariant
    

'''
for i in range (1, n_t):
    #find boundary condition for x = x_o from the derivative flux condition - i.e. find how the h(x_o) evloves in time
    h_xt[i][0] = h_xt[i-1][0] -  h_t*del_t

for i in range (0, n_t): 
    for j in range (1, n_x): #spatial index
        h_xt[i][j] = h_xt[i][j-1] + h_x[i]*del_x #evolves slope along the same curve

'''
'''
for i in range (0, n_t-1):
    for k in range (1, n_t-1):
        h_x_t[k-1] = (h_x[k] - h_x[k-1])/(del_t)
    for j in range (0, n_x):
        h_xt[i+1][j] =  h_xt[i][j] + del_t*(np.sign(h_x[i])*np.power(abs(h_x[i]),6/5) -2/5*np.sign(h_x[i])*np.power(abs(h_x[i]),-6/5)*h_x_t[i])
        #double check that this is true 
'''
fig, ax = plt.subplots()
plt_der1 = ax.plot(t, h_x[0], '--', label = "Numerical approximation with Forward Euler")
plt_der = ax.plot(t, h_x[1], '--', label = "Numerical approximation with Forward Euler")
plt_der2 = ax.plot(t, h_x[2], '--', label = "Numerical approximation with Forward Euler")
plt_der3 = ax.plot(t, h_x[3], '--', label = "Numerical approximation with Forward Euler")
#plt_i = ax.plot(x, h_xt[0], label = "initial profile")
#plt_xi = ax.plot(t, xi, label = "xi")
#plt_2 = ax.plot(x, h_xt[int(n_t/10)], label = "initial profile")
#plt_3 = ax.plot(x, h_xt[int(n_t/5)], label = "initial profile")
#plt_half = ax.plot(x, h_xt[int(n_t/2)], label = "profile2")
#plt_final = ax.plot(x, h_xt[int(n_t)-1], label = "initial profile")

#ax.legend( ["Numerical approximation of h_x"])
ax.legend( ["a = 0.01","a = 0.1", "a=1", "a = 2"] )

plt.xlabel("t [dimensionless]")
plt.ylabel("h_x [dimensionless]")
plt.title("Euler Fwd on h_x")
plt.show()