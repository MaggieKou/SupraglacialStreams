#adiabatic river case using forward euler method 
import numpy as np
from matplotlib import pyplot as plt 
import math

#numerical method to implement equation 21 from Moulins.pdf by Neufeld and Devauchelle

n = 10000 #num of steps 
#assuming t0 = 0
tf = 10
del_t = tf/n #step size, assuming uniform
H_o = 2 #initial H
H_prime_o = 5 #initial derivative

t = np.linspace(0, tf, n)
H = np.zeros([n])
H[0] = H_o
H_prime = np.zeros([n])
H_prime[0] = H_prime_o

def fwd_euler(t,h):
    for i in range (1, n):
        h[i] = h[i-1] + del_t*5/2*h**(11/5)*(1-h)
    return h

H_prime = fwd_euler(t, H_prime)

for i in range (1,n):
    H[i] = H[i-1] + del_t*H_prime[i]
print (t , "\n" , H)

fig, ax = plt.subplots()
plt_approx = ax.plot(x, H, '--', label = "Numerical approximation with Forward Euler")
ax.legend( ["Numerical approximation with Forward Euler Method"])

plt.xlabel("position [dimensionless]")
plt.ylabel("depth [dimensionless]")
plt.title("Euler Fwd on Receding River Case")
plt.show()


    
     
