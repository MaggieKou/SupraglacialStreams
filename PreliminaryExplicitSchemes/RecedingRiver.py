#receding river case using a forward euler method, as before 
import numpy as np
from matplotlib import pyplot as plt 
import math

#numerical method to implement from Moulins.pdf by Neufeld and Devauchelle

n = 10000000 #num of steps 
xo = -5000
xf = 1000
del_t = xf/n #step size, assuming uniform
H_o = 5 #initial H
H_prime_o = 0.4 #initial derivative

t = np.linspace(xo, xf, n)
H = np.zeros([n])
H[0] = H_o
H_prime = np.zeros([n])
H_prime[0] = H_prime_o

def fwd_euler(t,h):
    for i in range (1, n):
        h[i] = h[i-1] + del_t*5/2*h[i-1]**(11/5)*(1-h[i-1]) #evaluates a single step 
    return h

H_prime = fwd_euler(t, H_prime)

for i in range (1,n):
    H[i] = H[i-1] + del_t/2*(H_prime[i]+H_prime[i-1])
print (t , "\n" , H)

fig, ax = plt.subplots()
plt_derivative = ax.plot(t, H_prime, label = 'Numerical Approximation of Derivative')
plt_approx = ax.plot(t, H, '--', label = "Numerical approximation with Forward Euler")
ax.legend( ['Numerical Approximation of Derivative', "Numerical approximation with Forward Euler Method"])

plt.xlabel("position [dimensionless]")
plt.ylabel("depth [dimensionless]")
plt.title("Euler Fwd on Receding River Case")
plt.show()
