# -*- coding: utf-8 -*-
"""
All the solvers needed to iterate over a given domain for the receding assuption solution
"""
import numpy as np 
def xi(S):
    xi = (S)/np.power(1 + np.power(abs(S), 2), 1/2)
    return xi
#solves for xi based on current slope values

def h_t(xi_a, xi_b, del_t):
    return -(np.power(abs(xi_b), 6/5, dtype = float) + 2/5*np.power(abs(xi_b), -6/5, dtype = float)/del_t*(xi_b-xi_a))
#h_t based on current xi, uses previous xi value as a parameter

def ODE_solver(H_p,v):
    #takes initial input of slope at a given point and spits out the concavity from the reduced ode
    H_pp = 5/2*np.sign(H_p)*np.power(abs(H_p),11/5, dtype=float)*np.power(1 + np.power(abs(H_p),2, dtype=float)*np.power(abs(v), 10, dtype=float),4/10, dtype=float)*(np.sign(H_p, dtype=float)*np.power(abs(H_p), 1/5, dtype=float) + np.power(1 + np.power(abs(H_p),2, dtype=float)*np.power(abs(v), 10, dtype=float),1/2, dtype=float))
    return H_pp

def init_BC(h, S, t, v):
    del_t = (t[len(t)-1] - t[0])/len(t)
    S_o = BC_solver(S[0][0], t, v)
    
    for j in range (1, len(t)): #slope at given point
        S[0][j] = S_o[j]
        h[0][j] = h[0][j-1] + del_t*h_t(xi(S[0][j-1]),xi(S[0][j]),del_t)
    
    return h, S
        
def BC_solver(S_in, t, v):
    #given initial slope and v, calculate slope for the rest of the time 
    del_t = (t[len(t)-1] - t[0])/len(t)

    S_o = np.zeros([len(t)])
    S_o[0] = S_in #initial condition
    for j in range(1, len(t)):
        S_o[j] = S_o[j-1] - del_t*np.power(v, 12)*ODE_solver(np.power(v, -5, dtype = float)*S_o[j-1], v)
    return S_o

def iter_max_thresh(h, S, t, x, err):
    #1st, son point 1,1, solve for the current
    del_t = (t[len(t)-1] - t[0])/len(t)
    del_x = (x[len(x)-1] - x[0])/len(x)
    
    for j in range(1,len(t)):
        for i in range(1, len(x)):
            h[i][j] = h[i][j-1] + del_t*h_t(xi(S[i][j-1]),xi(S[i-1][j]),del_t) #assume the slope of the adjacent point
            #use this to calculate the local slope 
            S[i][j] = (h[i][j] - h[i-1][j])/del_x
            #calculate
            xi_temp = xi(S[i][j])
            xi_temp2 = xi(S[i][j-1])
            while(abs(xi_temp - xi_temp2)/xi_temp > err):
                #repeat
                xi_temp2 = xi_temp
                h[i][j] = h[i][j-1] + del_t*h_t(xi(S[i][j-1]),xi_temp2,del_t) #assume the slope of the adjacent point
                #use this to calculate the local slope 
                S[i][j] = (h[i][j] - h[i-1][j])/del_x
                xi_temp = xi(S[i][j])
    return h,S
#threshold for maximal iteration (error provided by the calculated xi deviation)
