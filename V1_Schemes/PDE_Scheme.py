# -*- coding: utf-8 -*-
"""
Created on Tue May 30 12:01:58 2023

scheme for marching through both space and time 
"""
import numpy as np 

def iterator(Sxlast, Stlast, hxlast, htlast, del_x, del_t, ep):
    h = 0
    S_t = 0
    h_t = 0
    S = Sxlast
    while (abs(Sxlast - (h-hxlast)/del_x + Stlast - (h-htlast)/del_t)>ep):
    #old scheme
        S_t = (S - Stlast)/del_t
        S_x = (S - Sxlast)/del_x #superfluous for now
        
        h_t = -2/5*np.power(abs(S),-6/5)*S_t - np.power(abs(S),6/5)/(np.sqrt(1+S**2))
        h = htlast + del_t*h_t
        S = (h-hxlast)/del_x
        print("Stlast", Stlast, "S", S)
    return (h, S)

def marching(S, h, t, x, ep):
    "given two m x n matrices that store both S and h at each point"
    #assume that the S and h are both populated for the [0][j] column and the [i][0] row
    n_t = len(S[0]) #assuming a rectangular domain
    n_x = len(S)
    del_t = (t[n_t-1]-t[0])/n_t
    del_x = (x[n_x-1] - x[0])/n_x
    for j in range (1,n_t):
        for i in range (1, n_x):
            h[i][j], S[i][j] = iterator(S[i-1][j], S[i][j-1], h[i-1][j], h[i][j-1],del_x, del_t,ep)
            print("space step:",i, "\n time step:", j, "\n h:", h[i][j])
    return (S, h)
# populates h and S

def incising (S, h, t, x, ep):
    n_t = len(S[0]) #assuming a rectangular domain
    n_x = len(S)
    del_t = (t[n_t-1]-t[0])/n_t
    del_x = (x[n_x-1] - x[0])/n_x
    
    for j in range (1,n_t):
        for i in range (1, n_x):
            h[i][j] = h[i][j-1] + del_t*(np.power(abs(S[i-1][j-1]), 6/5))/np.sqrt(1+ S[i-1][j-1]**2)
            S[i][j] = (h[i][j] - h[i-1][j])/del_x 
    return (S, h)