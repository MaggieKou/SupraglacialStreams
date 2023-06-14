#takes Linear Approximation of Incising River Case to solve Explicitly 
import numpy as np 
from matplotlib import pyplot as plt

#ssuming constant slope case
F = 1373 #W/m^2
rho = 0.001 #density of freshwaterkg/m^3
L = 3*10**5 #latent heat of the ice
h_t = F/(rho*L) #slope 
h_x = h_t**(5/6)

#time stepping linear 
t_f =10
t_0 = 0
n_t = 100
del_t = (t_f - t_0)/n_t

#spatial steps
x_0 = 0
x_f = 5
n_x = 100
del_x = (x_f - x_0)/n_x
#IC
h_top = 10
a = 2 #slope
b = 1
h_in = h_top -a*np.linspace(x_0,x_f, n_x) #linear initial distribution
#solutions 
t = np.linspace(t_0,t_f, n_t) 
x = np.linspace(x_0,x_f, n_x) 

h_xt =np.zeros([n_x,n_t])
h_xt[0] = h_in


for i in range (1, n_t): #timestep in t 
    for j in range (0, n_x):
        h_xt[i][j] = h_xt[i-1][j] - del_t*((a+b*x[j]**(6/5))) #let's come back to this scheme

fig, ax = plt.subplots()
#plt_approx = ax.plot(t, H_prime, '--', label = "Numerical approximation with Forward Euler")
plt_i = ax.plot(x, h_xt[0], label = "initial profile")
plt_half = ax.plot(x, h_xt[int(n_t/2)], label = "profile2")
plt_final = ax.plot(x, h_xt[int(n_t)-1], label = "initial profile")

#ax.legend( ["Numerical approximation of H'", "Numerical approximation of H"])
ax.legend( ["Initial Profile","Profile after some time, t/4","Profile after some time, t/2", "Profile after some time, t"])

plt.xlabel("x [dimensionless]")
plt.ylabel("h [dimensionless]")
plt.title("Euler Fwd on h")
plt.show()