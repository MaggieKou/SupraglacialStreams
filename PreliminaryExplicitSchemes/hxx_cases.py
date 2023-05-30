# -*- coding: utf-8 -*-
'''import numpy as np
import math
from matplotlib import pyplot as plt 


#ssuming constant slope case
F = 1373 #W/m^2
rho = 0.001 #density of freshwaterkg/m^3
L = 3*10**5 #latent heat of the ice
h_t = F/(rho*L) #slope 


#with scheme, solve after linearization
t_f = 20
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
So = 0.1 #slope
cases for hxx where hxtx varies as a function of hxt, assuming that hxx is nonzero
'''


h_in = h_top -So*np.linspace(x_0,x_f,n_x) #linear initial distribution
#solutions 
t = np.linspace(t_0,t_f, n_t) 
x = np.linspace(x_0,x_f, n_x) 


hxx = np.zeros([n_t,n_x]) #Potentially varies in both time and space
hx = np.zeros([n_t,n_x])
hxt = np.zeros([n_t,n_x])
a = np.zeros([n_t])

h = np.zeros([n_t,n_x])
h[0] = h_in

'''
def linear (hx, hxt, a): # where hxt is in the form a(t)x, then hxtx = a(t)
    for i in range (0, n_t):
        for j in range (0, n_x):
            hxx[i][j] = a[i]/(hx[i][j]*(3*hx[i][j]**(2/5) + 6/5*hxt[i][j]))
    return hxx
'''
#time step for hx
xi = np.zeros([n_t])
xi[0] = -np.power(So,-1/5) 

for i in range (1, n_t):
    xi[i]= xi[i-1] + del_t*(h_t/5 -1/2*(xi[i-1])**(-6)) #gives xi which is h_x^(-1/5)

for i in range (0, n_t):
    hx[i][0] = xi[i] #initialized time steps for h_x at initial position 

#initialize hxt
for i in range (0, n_t): #initial derivative
    hxt[i][0] = (5/(2*xi[i]) - h_t)*np.power(hx[i][0],6/5)


for i in range (0, n_t):
    for j in range (1, n_x):
        hx[i][j] = hx[i][j-1] +a[i]*del_x/(hx[i][j-1]*(3*np.sign(hx[i][j-1])*np.power(hx[i][j-1],2/5) + 6/5*hxt[i][j-1])) #no operation on time, forward stepping the position

fig, ax = plt.subplots()

plt_der1 = ax.plot(x, hx[0], '--', label = "Numerical approximation with Forward Euler")
plt_der2 = ax.plot(x, hx[int(n_t/2)], '--', label = "Numerical approximation with Forward Euler")
plt_der3 = ax.plot(x, hx[n_t -1], '--', label = "Numerical approximation with Forward Euler")

#time and space step for hx
ax.legend( ["initial","tf/2", "tf"])
#ax.legend( ["a = 0.1","a = 1", "a=2", "a = 3"] )

plt.xlabel("t [dimensionless]")
plt.ylabel("h_x [dimensionless]")
plt.title("Euler Fwd on h_x")
plt.show()
