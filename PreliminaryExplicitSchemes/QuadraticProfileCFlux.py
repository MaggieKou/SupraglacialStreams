''' Solves the linear profile for the quadratic case
'''
import numpy as np
import math
from matplotlib import pyplot as plt 


#ssuming constant slope case
F = 1373 #W/m^2
rho = 0.001 #density of freshwaterkg/m^3
L = 3*10**5 #latent heat of the ice
h_t = -F/(rho*L) #slope 


#with scheme, solve after linearization
t_f = 2
t_0 = 0
n_t = 100
del_t = (t_f - t_0)/n_t

#spatial steps
x_0 = 0.1
x_f = 2
n_x = 10
del_x = (x_f - x_0)/n_x


t = np.linspace(t_0,t_f, n_t) 
x = np.linspace(x_0,x_f, n_x) 


#initial conditions - assuming a quadratic profile of the form h = h_top - C*x^2 + Skew*x
h_top = 500
C = 4 #concavity 
Skew = 0.5

h_o = h_top - C*x**2 + Skew*x

h =np.zeros([n_t,n_x]) #solution space
h[0] = h_o

#time step for hx
xi = np.zeros([n_t,n_x])
hx = np.zeros([n_t,n_x])
hxt = np.zeros([n_t,n_x])

for i in range (0, n_x):
    xi[0][i] = -np.sign(2*C*x[i] - Skew)*np.power(abs(2*C*x[i] - Skew),-1/5) #initial u sub, dependent on spatial derivative
    hx[0][i] = np.power(xi[0][i],-5)
    
for i in range (1, n_t):
    for j in range (0, n_x):
        xi[i][j]= xi[i-1][j] + del_t*(h_t/2 +1/2*(xi[i-1][j])**(-6)) #gives xi which is h_x^(-1/5)
        hx[i][j] = np.power((xi[i][j]),-5)
        
for i in range (1, n_t):
    #find boundary condition for x = x_o from the derivative flux condition - i.e. find how the h(x_o) evloves in time
    h[i][0] = h[i-1][0] - h_t*del_t

for i in range (1, n_t):
    for k in range (0, n_x):
        hxt[i][k] = (hx[i][k] - hx[i-1][k])/(del_t) #gives slope
    
for i in range (1, n_t):
    for j in range (1, n_x): #spatial index
    #for all time, evolve from the initial height at xo for that time!
        h[i][j] =  h[i-1][j] + del_t*(-np.power(abs(hx[i-1][j]),6/5) +2/5*np.power(abs(hx[i-1][j]),-6/5)*hxt[i][j])
        
fig, ax = plt.subplots()
plt_i = ax.plot(x, hxt[0], label = "initial profile")
plt_i = ax.plot(x, hx[0], label = "initial profile")

plt_i = ax.plot(x, xi[0], label = "initial profile")
plt_i = ax.plot(x, xi[1], label = "initial profile")
#plt_i = ax.plot(x, xi[99], label = "initial profile")
#plt_i = ax.plot(x, xi[100], label = "initial profile")


'''
plt_i = ax.plot(x, xi[0], label = "initial profile")
plt_i = ax.plot(x, xi[1], label = "initial profile")
plt_i = ax.plot(x, xi[2], label = "initial profile")
#plt_i = ax.plot(x, xi[100], label = "initial profile")

plt_i = ax.plot(x, hx[0], label = "initial profile")
plt_i = ax.plot(x, hx[1], label = "initial profile")
#plt_i = ax.plot(x, hx[10], label = "initial profile")
plt_i = ax.plot(x, hx[int(t_f/2)], label = "initial profile")
plt_i = ax.plot(x, hx[int(t_f)-1], label = "initial profile")

plt_i = ax.plot(x, h[0], label = "initial profile")
#plt_i = ax.plot(x, h[100], label = "initial profile")
plt_i = ax.plot(x, h[1], label = "initial profile")
#plt_i = ax.plot(x, h[10], label = "initial profile") # fix code, not sure what's happening here
plt_i = ax.plot(x, h[int(t_f)-1], label = "initial profile")


plt_i = ax.plot(x, h[0], label = "initial profile")
#plt_i = ax.plot(x, h[100], label = "initial profile")
plt_i = ax.plot(x, h[1], label = "initial profile")
#plt_i = ax.plot(x, h[10], label = "initial profile") # fix code, not sure what's happening here
#plt_i = ax.plot(x, h[int(t_f)-1], label = "initial profile")
'''
ax.legend( ["h(initial)", "h shortly after", "h after some dimensionless time", "h stable"])
#ax.legend( ["a = 0.1","a = 1", "a=2", "a = 3"] )

plt.xlabel("x [dimensionless]")
plt.ylabel("h [dimensionless]")
plt.title("Euler Fwd on h for the quadratic case")
plt.show()