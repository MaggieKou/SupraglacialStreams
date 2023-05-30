# -*- coding: utf-8 -*-
"""
Created on Tue May 30 10:07:47 2023
Solves the boundary condition ODE with an explicit forward euler scheme: \frac{7 S_{t}^2}{5S} - S_{tt} = 0 |_{x = x_o}
"""
import numpy as np

"constants that set initial conditions"
F = 1373 #W/m^2
rho = 0.001 #density of freshwaterkg/m^3
L = 3*10**5 #latent heat of the ice
Q = F/(rho*L) #slope 
g = 9.81 
Cf = 0.04


#iteratively solves for the time initial slope at a fixed point
#since this technically applies for all x, wee are reducing our effective dimensions and can apply this in time to all of the initial conditions to solve
def initialize(h, S, t, x, ep, R): #R is the aspect ratio
    n_t = len(S[0])
    del_t = (t[n_t-1]-t[0])/n_t
    #sets the solution output space
    St = np.zeros([n_t]) 
    
    So = S[0][0] #from the initial conditions
    
    #sets initial time derivative of the slope 
    Sto = -5/2*(Q**(1/5)*g**(2/5))/(R**(1/5)*Cf**(2/5)*np.sign(So)*np.power(abs(So), -7/5))
    St[0] = Sto
    
    #solves for all of S(x_o)
    S[0] = solve_S(n_t, S[0], St, del_t, ep)
    h[0] = solve_h(h[0], S[0],del_t, n_t, ep)
    return (S,h)
    
    
def solve_S(n_t, S, St, del_t, ep):
    Snw = 0
    for i in range (1, n_t):
        St[i] = St[i-1] + 7/5*del_t*(St[i-1]**2)/S[i-1] #solves for the derivative at this stage
        S[i] = S[i-1] + del_t*(St[i]+St[i-1])/2
        while (abs(Snw-St[i]) > ep):
            Snw = St[i-1] + 7/5*del_t*(St[i-1]**2)/S[i-1] #recalculates 
            S[i] = S[i-1] + del_t*(Snw+St[i-1])/2
        St[i] = Snw
    return S

def solve_h(h, S, del_t, n_t, ep):
    for i in range (1, n_t): 
        h[i] = h[i-1] + del_t*(-2/5*np.power(abs(S[i]),-6/5))*((S[i]-S[i-1])/del_t) - (np.power(abs(S[i]), 6/5))*(np.sqrt(1+S[i]**2))
    return h