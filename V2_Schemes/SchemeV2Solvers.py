# -*- coding: utf-8 -*-
"""
Solvers for V2 of the explicit scheme.
IC: the initial topography of the ice. 
BC: the width as a function of time for the first point as a boundary condition. This is a prescribed melting pattern at some point.  
"""
import numpy as np 
from scipy.signal import savgol_filter
import math

"helper functions"

def solve_St(S_vec, del_t, i): #given some S vector i.e. nx1 list for time, and some position, i 
    if(i>1 and i <len(S_vec) -1):
        St = (S_vec[i+1] - S_vec[i-1])/(2*del_t)
    elif i == 1: 
        St = (S_vec[1]-S_vec[0])/del_t
    else: 
        St = (S_vec[i] - S_vec[i-1])/del_t
        
        #leaves the initial St value empty 
    return St

def solve_ht(S, St):
    #IN: S value and St value (guessed)
    #OUT: ht value
    #this is a scalar function 
    ht = -2/5*np.power(abs(S), -6/5)*St - np.power(abs(S), 6/5)/(np.sqrt(1 + abs(S)**2))
    return ht



'''functions all for the first spatial point'''
def S_BC(S, w_t, del_t, Q, R):
    #given BC for the width and the solar heating, and the initial IC slope, the rate of change of the width as a function of time, 
    #calculates the slope at that point
    g = 9.81
    Cf = 0.1 
    
    w0 = -np.sign(g*S[0][0]/(Cf*R*Q**2))*np.power(np.abs(g*S[0][0]/(Cf*R*Q**2)), -1/5) 
    w = np.zeros([len(S[0])])
    w[0] = w0
    
    for j in range (1, len(w)):
        w[j] = w[j-1] + del_t*w_t[j] #initializes width with a taylor expnasion 
        S[0][j] = -Cf*R*(Q**2)/(g*np.power(abs(w[j]),5)) #shouldn't have any sign errors as the w >0 always 
        #print ("w0", w0, "w_t", w_t[j], "w", w[j], "S", S[0][j])
    return S


def h_BC(del_t, h, S):
    #IN: S for IC and BC, h for IC
    #OUT: h for BC
    for j in range (1, len(h[0])):
        h[0][j] = h[0][j-1] + del_t*solve_ht(S[0][j], solve_St(S[0], del_t, j)) #use information at this point
    return h

"functions for stepping from the ICs and BCs"
def solve_fountain(h, S, i, j, del_t, del_x, ep, Q):
    eta = 0.01
    L = 3.35*10**(5)
    g = 9.8 
    Const = 1/2*(math.pi/(2*eta))**(3/8)*(g/L)
    ht_temp = Const*np.power(abs(S[i-1][j]), 19/16)*np.power(abs(Q), 5/8)
    h[i][j] = h[i][j-1] + del_t*ht_temp #guess h of new point
    S_temp = (h[i][j]- h[i-1][j])/del_x #compare with slope 
    S[i][j] = S[i-1][j]
    
    while (abs(S[i][j] - S_temp)/S_temp> ep):
        S[i][j] = S_temp
        ht_temp = Const*np.power(abs(S[i][j]), 19/16)*np.power(abs(Q), 5/8)
        h[i][j] = h[i][j-1] + del_t*ht_temp
        S_temp = (h[i][j]- h[i-1][j])/del_x
        print("S_temp", S_temp, "S", S[i][j],"err", abs(S_temp - S[i][j])/S_temp)

    S[i][j] = S_temp
    return h[i][j], S[i][j]

    
def solve_incising(h, S, i, j, del_t, del_x, ep):
    #IN: slope matrix, use S from previous spatial step to guess the LHS, then iterate by guessing h frome PDE, find slope
    ht_temp = -np.power(abs(S[i-1][j]),6/5)/(np.sqrt(1+np.power(abs(S[i-1][j]),2))) #guess based on previous timestep slope
    h[i][j] = h[i][j-1] + del_t*ht_temp #guess h of new point
    S_temp = (h[i][j]- h[i-1][j])/del_x #compare with slope 
    S[i][j] = S[i-1][j]
    
    while (abs(S[i][j] - S_temp)/S_temp> ep):
        S[i][j] = S_temp
        ht_temp = -np.power(abs(S[i][j]),6/5)/(np.sqrt(1+np.power(abs(S[i][j]),2)))
        h[i][j] = h[i][j-1] + del_t*ht_temp
        S_temp = (h[i][j]- h[i-1][j])/del_x
        print("S_temp", S_temp, "S", S[i][j],"err", abs(S_temp - S[i][j])/S_temp)

    S[i][j] = S_temp
    return h[i][j], S[i][j]

def solve_widening(h, S, i, j, del_t, del_x, ep):
    #IN: slope matrix, use S from previous spatial step to guess what the new slope is, 
    #then use slope to guess point from previous space step, iterate based on slope - this convergence is not guaranteed
    #assume that spatial slope is more smooth than temporal
    S[i][j] = S[i][j-1] - del_t*np.sign(S[i-1][j])*np.power(abs(S[i-1][j]),7/5)/(np.sqrt(1 + np.power(abs(S[i-1][j]))))
    h[i][j] = h[i-1][j] + del_x*S[i][j]
    return h[i][j], S[i][j] 

def analytic_adiabatic(h, S, i, j, x, del_x):
    #provides analytic solution to the adiabatic case 
    h[i][j] = h[0][j] - 2*np.power(3, 1/6)*np.power(abs(x[0] - x[i]), 1/6)
    S[i][j] = (h[i][j]- h[i-1][j])/del_x
    print("aadiabatic solutions:", "h", h[i][j], "j", j)
    return h[i][j], S[i][j]

def strang_split(h, S, n_x, n_t, del_t, del_x, x, err):
    #IN: h and S matrices populated with BCs and ICs. 
    h_1 = h
    S_1 = S
    
    #first, cacluate the two separate h points for one time and space step 
    for j in range (1, n_t):
        for i in range (1, n_x):
            h_1[i][j], S_1[i][j] = analytic_adiabatic(h, S, i, j, x, del_x) #first, approximate point with adiabatic case
            h[i][j], S[i][j] = solve_incising(h_1, S_1, i, j, del_t, del_x, err) #uses estimate from adiabatic case
    
    return solve_spacetime(h, S, n_x, n_t, del_t, del_x, err)
def solve_spacetime(h, S, n_x, n_t, del_t, del_x, err):
    #IN: mxn matrix for h with ICs and BCs populated, mxn S, error threshold
    #OUT: mxn matrix, solved
    for j in range(1, n_t):
        for i in range(1, n_x):
            #start with the first point (1, 1)
            h[i][j] = h[i][j-1] + del_t*solve_ht(S[i][j], solve_St(S[i], del_t, i)) #guess that we have the slope as the previous point
            S[i][j] = (h[i][j] - h[i-1][j])/del_x #linearize the slope
            print("first guess", S[i][j])
            h_t_temp = solve_ht(S[i][j], (S[i][j] - S[i][j-1])/del_t) #use current point to calculate h_t
            #h[i][j] = h[i][j-1] + del_t*h_t_temp
            print("h_t_temp", h_t_temp, "h", h[i][j],"S", S[i][j], "err", err)
            #h[i][j], S[i][j] = error_iterator(h_t_temp, err, h, S, i, j, del_t, del_x)
    return h
            
    #using the S and h values from both time and space, attempt to converge on an S value 
    

def smoothing (h, S):
    S = savgol_filter(S,1000,4,axis=0)
    h = savgol_filter(h,300,6,axis=0)
    return h,S

def error_iterator(h_t_temp, err, h, S, i, j, del_t, del_x): 
    h2 = (h[i][j] - h[i][j-1])/del_t
    while(abs(h2 - h_t_temp)< err):
        h2 = h_t_temp
        h[i][j] = h[i][j-1] + del_t*h_t_temp
        S[i][j] = (h[i][j] - h[i-1][j])/del_x #reclaculates slope
        h_t_temp = solve_ht(S[i][j], solve_St(S[i], del_t, i)) #use current point to calculate h_t again
    return h[i][j], S[i][j]