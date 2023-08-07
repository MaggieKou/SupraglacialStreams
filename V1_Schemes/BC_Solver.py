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

def corrector(h, S, j, del_t, del_x, ep):
    #corrects based on a fully filled row 
    n_x = len(h)
    for i in range (1, n_x-1):
        S_temp = (h[i+1][j]-h[i-1][j])/(2*del_x)
        while (abs(S_temp - S[i][j]) > ep):
            S[i][j] = S_temp
            h_t = -(np.power(abs(S[i][j]),6/5))/(np.sqrt(1 + (S[i][j])**2)) -2/5*np.power(abs(S[i][j]), -6/5)*(S[i][j]-S[i][j-1])/del_t
            h[i][j] = h[i][j-1] + del_t*h_t
            S_temp = (h[i+1][j]-h[i-1][j])/(2*del_x)
        S[i][j] = S_temp
    S[i][j] = (h[i+1][j]-h[i-1][j])/(2*del_x)
    return h, S
        
def iter_row(h, S, j, start_i, end_i, del_t, del_x, ep): # j >0
    S_1 = S[start_i][int(j)-1]#guess 
    h_t = -(np.power(abs(S_1),6/5))/(np.sqrt(1 + (S_1)**2)) #guess no time change in slope for now 
    h[int(start_i)][int(j)] = h[int(start_i)][int(j)-1] + del_t*h_t #takes the same slope as the time point before it
    
    for i in range (start_i + 1,end_i):  
        print(i,j)
        if (i > 2): #calculates based on adjacent slope
            h_t = -(np.power(abs(S[i-1][j]),6/5))/(np.sqrt(1 + np.power(abs(S[i][j]),2))) -2/5*np.power(abs(S[i-1][j]), -6/5)*(S[i-1][j]-S[i][j-1])/del_t

        else: #if not first or second point 
            h_t = -(np.power(abs(S[i][int(j-1)]),6/5))/(np.sqrt(1 + np.power(abs(S[i][j]),2))) #just guess no time dependence for now, we'll use the slope from the previous point to guess
        
        h[i][int(j)] = h[i][int(j-1)] + del_t*h_t
        S[i][int(j)] = (h[i][int(j)]-h[i-1][int(j)])/del_x #this is the slope using the point next to it
        
        #recalc
        h_t = -(np.power(abs(S[i][j]),6/5))/(np.sqrt(1 + np.power(abs(S[i][j]),2))) -2/5*np.power(abs(S[i][j]), -6/5)*(S[i][j]-S[i][j-1])/del_t
        h[i][j] = h[i][j-1] + del_t*h_t #realculates based on previous point 
        
        S_temp = (h[i][j]-h[i-1][j])/del_x #recalc slope 
        
        while (abs(S_temp - S[i][j])>ep):
            S[i][j] = S_temp
            h_t = -(np.power(abs(S[i][j]),6/5))/(np.sqrt(1 + (S[i][j])**2)) -2/5*np.power(abs(S[i][j]), -6/5)*(S[i][j]-S[i][j-1])/del_t
            h[i][j] = h[i][j-1] + del_t*h_t
            S_temp = (h[i][j]-h[i-1][j])/del_x
        S[i][j] = S_temp
    S[start_i][j] = (h[start_i +1][j]-h[start_i][j])/del_x
    
    return corrector(h, S, j, del_t, del_x, ep)
        

#iteratively solves for the time initial slope at a fixed point
#since this technically applies for all x, wee are reducing our effective dimensions and can apply this in time to all of the initial conditions to solve
def initialize(h, S, t, x, ep, R, w): #R is the aspect ratio
    n_t = len(S[0]) #BC conditions applied to this column
    n_x = len(S)
    del_x = (x[n_x-1]-t[0])/n_x
    del_t = (t[n_t-1]-t[0])/n_t
    #sets the solution output space
    St = np.zeros([n_t]) 
    
    #we want to compute the slope at this point, directly proporitional to the initial w
    #given h[0][0] from initial conditions, and S[0][0], the slope from the neighbouring IC
    So = (Cf*R*Q**2)/(g*w**5) #this is our initial slope
    for j in range (1, n_t):
        S[0][1] = So
    
    for j in range (2, n_t):
        h,S = iter_row(h, S,j, 1 , n_x, del_t, del_x, ep)
         
    # we will iterate on these once the other x's have been filled out too
    
    #R is prescribed by user
    return (S,h)
    
    
def solve_S(n_t, S, St, del_t, ep):
    Snw = 0 #new slope time derivative, to compare against
    for i in range (1, n_t):
        St[i] = St[i-1] + 7/5*del_t*((St[i-1])**2)/S[i-1] #solves for the derivative at this stage from previous "S" value
        #S[i] = S[i-1] + del_t*(St[i]+St[i-1])/2 #averages the left and right slopes 
        S[i] = S[i-1] + del_t*St[i]#new slope at this poin 
        while (abs(Snw-St[i]) > ep):
            St[i] = Snw
            Snw = St[i-1] + 7/5*del_t*(St[i]**2)/S[i] #recalculates with new S
            S[i] = S[i-1] + del_t*Snw #avgs
            print(St[i])
        St[i] = Snw
    return S

def solve_h(h, S, del_t, n_t, ep):
    for i in range (1, n_t): 
        h[i] = h[i-1] + del_t*(-np.power(abs(S[i]),6/5)/np.sqrt(1+S[i]**2))
        #h[i] = h[i-1] + del_t*(-2/5*np.power(abs(S[i]),-6/5))*((S[i]-S[i-1])/del_t) - (np.power(abs(S[i]), 6/5))*(np.sqrt(1+S[i]**2))
    return h