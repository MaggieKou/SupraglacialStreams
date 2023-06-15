# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np

def ODE_solver(H_p,v):
    #takes initial input of slope at a given point and spits out the concavity from the reduced ode
    H_pp = -5/2*np.sign(H_p, dtype=float)*np.power(abs(H_p),11/5, dtype=float)*(np.sign(H_p, dtype=float)*np.power(abs(H_p), 1/5, dtype=float)/(np.sqrt(1 + np.power(abs(H_p),2, dtype=float)*np.power(abs(v), 10, dtype=float))) + 1)
    return H_pp
def timestepper(H_in, x, t, v):
    del_x = (x[len(x)-1]-x[0])/len(x)
    del_t = (t[len(t)-1]-t[0])/len(t)
    #assume that slope is relatively continuous to extract for slopes near boundaries
    H = np.zeros([len(x),len(t)]) #index space first, then index time
    S = np.zeros([len(x),len(t)]) #index space first, then index time
    
    #initialize interior points
    for i in range (1, len(x)):
        H[i][0] = H_in[i]
        S[i][0] = (H[i][0] - H[i-1][0])/del_x
    
    #initialize boundary points
    H[0][0] = H_in[0]
    S[0][0] = S[1][0]
    
    #populate based on receding slope
    for j in range(0, len(t)):
        for i in range (1, len(x)):
            H_p = S[i][j-1] + del_x*ODE_solver(S[i][j-1],v)
            H[i][j] = H[i][j-1] + del_t*H_p*np.power(abs(v), 6, dtype=float)
            if i == 1:
                S[0][j] = (H[1][j] - H[0][j])/del_x
            if i > 0:
                S[i][j] = (H[i][j] - H[i-1][j])/del_x
    return H, S