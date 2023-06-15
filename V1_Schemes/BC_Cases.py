# -*- coding: utf-8 -*-
"""
Created on Fri Jun  2 12:27:16 2023

@author: maggi
"""
import numpy as np

def S_x0_solver(t, S, w, R):
    #given some function w(t), this gives S(x0, t) for width at the initial point
    F = 1373 #W/m^2
    rho = 0.001 #density of freshwaterkg/m^3
    L = 3*10**5 #latent heat of the ice
    Q = F/(rho*L) #slope 
    g = 9.81 
    Cf = 0.04
    w_o = np.power(abs(Cf*R**3*Q**2/g), 1/5)
    n_t = len(t)
    del_t = (t[n_t-1] - t[0])/n_t
    for j in range(1, len(w)):
        S[0][j] = S[0][j-1] -5*S[0][j-1]*(w[j]-w[j-1])/w[j]
    return S

def setup(x, t, S, h, R, w, ep):
    #assume that initial conditions are already loaded in, with slope at initial point, assume one-sided continuity
    n_x = len(S)
    n_t = len(S[0])
    del_x = (x[n_x-1] - x[0])/n_x
    del_t = (t[n_x-1] - t[0])/n_t
    "1) Solve for BC slope at some initial point"
    S = S_x0_solver(t, S, w, R)
    "2) Solve for the next time step for BC based on the PDE"
    for j in range(1, n_t):
        h_t = -(np.power(abs(S[0][j]),6/5))/(np.sqrt(1 + (S[0][j])**2)) -2/5*np.power(abs(S[0][j]), -6/5)*(S[0][j]-S[0][j-1])/del_t
        h[0][j] = h[0][j-1] + h_t #temporary
    
    #initialized h and S points along finite perimeter
    return step_uniform(x,t,n_x, n_t, del_x, del_t, S, h, ep, w, R)

def step_uniform(x,t,n_x, n_t, del_x, del_t, S, h, ep, w, R):
    #given some function w(t), this gives S(x0, t) for width at the initial point
    F = 1373 #W/m^2
    rho = 0.001 #density of freshwaterkg/m^3
    L = 3*10**5 #latent heat of the ice
    Q = F/(rho*L) #slope 
    g = 9.81 
    Cf = 0.04
    w_o = np.power(abs(Cf*R**3*Q**2/g), 1/5)
    n_t = len(t)
    del_t = (t[n_t-1] - t[0])/n_t
    for j in range (1, n_t):
        for i in range(1, n_x):
             S[i][j] = S[i][j-1] -5*S[i][j-1]*(w[j]-w[j-1])/w[j]
             h_t = -(np.power(abs(S[i-1][j]),6/5))/(np.sqrt(1 + (S[i-1][j])**2)) -2/5*np.power(abs(S[i-1][j]), -6/5)*(S[i-1][j]-S[i][j-1])/del_t
             h[i][j] = h[i][j-1] + h_t*del_t
             #implement correcting scheme here
             
    return correction_scheme(x, t, n_x, n_t, del_x, del_t, S, h, ep)
             
def correction_scheme (x, t, n_x, n_t, del_x, del_t, S, h, ep):
    for i in range (1, n_x-2):
        for j in range (1, n_t):
           
            S_temp = (h[i+1][j] - h[i-1][j])/(2*del_x)
            while (abs(S[i][j]-S_temp)> ep):
                S[i][j] = S_temp
                h_t = -(np.power(abs(S[i][j]),6/5))/(np.sqrt(1 + (S[i][j])**2)) -2/5*np.power(abs(S[i][j]), -6/5)*(S[i][j]-S[i][j-1])/del_t
                h[i][j] = h[i][j-1] + h_t*del_t
                S_temp = (h[i][j] - h[i-1][j])/del_x
            S[i][j] = S_temp
    return h, S
        
def iterate_over(x,t,n_x, n_t, del_x, del_t, S ,h, ep):
    for j in range (1, n_t):
        for i in range (1, n_x):
            #initial guess for slope is the one before, in space
            h_t = -(np.power(abs(S[i-1][j]),6/5))/(np.sqrt(1 + (S[i-1][j])**2)) -2/5*np.power(abs(S[i-1][j]), -6/5)*(S[i-1][j]-S[i][j-1])/del_t
            h[i][j] = h[i][j-1] + h_t*del_t
           
            S_temp = (h[i][j] - h[i-1][j])/del_x
            while (abs(S[i][j]-S_temp)> ep):
                S[i][j] = S_temp
                h_t = -(np.power(abs(S[i][j]),6/5))/(np.sqrt(1 + (S[i][j])**2)) -2/5*np.power(abs(S[i][j]), -6/5)*(S[i][j]-S[i][j-1])/del_t
                h[i][j] = h[i][j-1] + h_t*del_t
                S_temp = (h[i][j] - h[i-1][j])/del_x
            S[i][j] = S_temp
    return h, S