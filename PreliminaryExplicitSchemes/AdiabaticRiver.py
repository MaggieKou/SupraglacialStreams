#adiabatic river case using forward euler method 
import numpy as np
from matplotlib import pyplot as plt 
import math

#numerical method to implement equation 21 from Moulins.pdf by Neufeld and Devauchelle

n = 100000 #num of steps 
x_o = 1 #initial position
xf = 100
del_x = (xf-x_o)/n #step size, assuming uniform
h_o = 100 #initial slope 


x = np.linspace(0, xf, n)
h = np.zeros([n])
h[0] =h_o
ep = 10**(-6) #perturbing factor
#plot logarthmic 
def fwd_euler(x,h):
    for i in range (1, n):
        h[i] = h[i-1] + del_x*((2**5)/((h_o - h[i-1]+ep)**5))
    return h

#plot analytic solution 
h_analytic = np.zeros([n])
for i in range (0,n):
    h_analytic[i] = h_o -2*3**(1/6)*(-x[i]+x_o)**(1/6)

h_approx = fwd_euler(x, h)
print (x , "\n" , h)

fig, ax = plt.subplots()
plt_analytic = ax.plot(x, h_analytic, '-', label="Analytic Solution")
plt_approx = ax.plot(x, h_approx, '--', label = "Numerical approximation with Forward Euler")
ax.legend( ["Analytic Solution","Numerical approximation with Forward Euler Method"])

plt.xlabel("position [dimensionless]")
plt.ylabel("depth [dimensionless]")
plt.title("Euler Fwd on Adiabatic River Case")
plt.show()

    
     
