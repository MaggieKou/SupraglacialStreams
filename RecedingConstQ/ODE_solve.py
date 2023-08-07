# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np

def ODE_solver(H_p,v):
    #takes initial input of slope at a given point and spits out the concavity from the reduced ode
    H_pp = -5/2*np.sign(H_p, dtype=float)*np.power(abs(H_p),11/5, dtype=float)*np.power(1 + np.power(abs(H_p),2, dtype=float)*np.power(abs(v), 10, dtype=float),3/10, dtype=float)*(np.sign(H_p, dtype=float)*np.power(abs(H_p), 1/5, dtype=float) + np.power(1 + np.power(abs(H_p),2, dtype=float)*np.power(abs(v), 10, dtype=float),6/5, dtype=float))
    return H_pp

def solve_row (H, S, del_x, v, t, del_t, i):
    #solves the solition for the time evolution at a specific position
    for j in range(1, len(t)):
            H_p = S[i][j-1] + del_t*np.power(v, 12)*ODE_solver(np.power(v, -5, dtype = float)*S[i][j-1],v) #steps the slope forwards based on what th slope was an instant ago
            H[i][j] = H[i][j-1] + del_t*H_p*np.power(v, 5)
            #solves for the rest of the row
            S[i][j] = H_p
    return H, S

def solve_slope(H, S, del_x, v, t, del_t, i):
    for j in range(1, len(t)):
        S[i][j]= S[i][j-1] + del_t*np.power(v, 12)*ODE_solver(np.power(v, -5, dtype = float)*S[i][j-1],v) #steps the slope forwards based on what th slope was an instant ago
            
    return S

def tranfr(h_in, x, t, v):
    del_t = (t[len(t)-1]-t[0])/len(t)
    #assume that slope is relatively continuous to extract for slopes near boundaries
    H = np.zeros([len(x),len(t)]) #index space first, then index time
    S = np.zeros([len(x),len(t)]) #index space first, then index time
    
    #scale the x-axis by velocity 
    x_new = np.zeros(len(x))
    for i in range(len(x)):
        x_new[i] = x[i]*np.power(abs(v), 6)
    
    del_x = (x_new[len(x)-1]-x_new[0])/len(x)

    #now plot the H value against the modified h value 
    for i in range(1, (x)):
        H[i][0] = v*h_in[i]  #sets the IC
        S[i][0] = np.power(v, -5, dtype = float)*(h_in[i] - h_in[i-1])/del_x
    #assuming that h_in is in the fixed frame
    H[0][0] = h_in[0]
    S[0][0] = S[1][0] 
    #print(S[0][0])
    
    "PLEASE FIX THIS PART"
    for i in range (0, len(x)):
        for j in range(1, len(t)):
            H_p = S[i][j-1] + del_t*np.power(v, 12)*ODE_solver(np.power(v, -5, dtype = float)*S[i][j-1],v) #steps the slope forwards based on what th slope was an instant ago
            H[i][j] = H[i][j-1] + del_t*H_p*np.power(v, 5)
            #solves for the rest of the row
            S[i][j] = H_p
    pass 
    
    
def timestepper(H_in, x, t, v): #assuming some initial profile, this outputs the dynamic solution in a fixed frame
    del_x = (x[len(x)-1]-x[0])/len(x)
    del_t = (t[len(t)-1]-t[0])/len(t)
    #assume that slope is relatively continuous to extract for slopes near boundaries
    H = np.zeros([len(x),len(t)]) #index space first, then index time
    S = np.zeros([len(x),len(t)]) #index space first, then index time
    
    #initialize interior points
    for i in range (1, len(x)):
        H[i][0] = H_in[i]
        S[i][0] = (H_in[i] - H_in[i-1])/del_x
    
    #initialize boundary points
    H[0][0] = H_in[0]
    S[0][0] = S[1][0] 
    #print(S[0][0])
    
    #populate based on receding slope
    for i in range (0, len(x)):
        H, S = solve_row(H, S, del_x, v, t, del_t, i)

    return H, S
#why is the initial slope so large?
'''
    for j in range(1, len(t)):
        for i in range (0, len(x)):
            H_p = S[i][j-1] + del_x*ODE_solver(S[i][j-1],v)
            H[i][j] = H[i][j-1] + del_t*H_p*np.power(abs(v), 6, dtype=float)
            if i == 1:
                S[0][j] = (H[1][j] - H[0][j])/del_x
            if i > 0:
                S[i][j] = (H[i][j] - H[i-1][j])/del_x
    return H, S
'''
