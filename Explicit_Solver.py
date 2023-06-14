#this python file is supposed to provide numerical solutions to a given functional input 
#this is for explicit schemes, in the form of dy/dx = f(x,y) for now
import numpy as np

def explicit_ode_solve(ode, x, BC):# BC is taken as y[0]
    #x is a n by 1 vector 
    n = len(x) 
    del_x = (x[n]-x[0])/n
    # y is a n by 1 vector 
    y = np.zeros([n])
    # the ode is an expression that relates x and y to its derivative, also a first, we call it as a function
    f = ode(x,y) #this will also return an n by 1 vector 
    y[0] = BC
    
    for i in range (1, n):
        y[n] = y[n-1] + del_x*f[n-1]
    return y
    
    