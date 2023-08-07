# -*- coding: utf-8 -*-
"""
Created on Tue May 30 12:25:52 2023

This file defines important global parameters for the problem setup
- currently broken: please fix this
"""
import numpy as np
import BC_Solver
import PDE_Scheme
import Plotter
import BC_Cases
from matplotlib import pyplot as plt
"Set the solver parameters here"

#TIMESTEPS
to = 0.1 #starting at the intial condition
tf = 1
n_t = 1000 #number of timesteps 
del_t = (tf-to)/n_t
#ho = np.linspace(to, tf, n_t) #initiate initial conditions 
t = np.linspace(to, tf, n_t)

#SPATIAL STEPS
xo = 0 #starting at the intial condition
xf = 1
n_x = 1000 #number of spatial steps
del_x = (xf-xo)/n_x
#ho = np.linspace(to, tf, n_t) #initiate initial conditions 
x = np.linspace(xo, xf, n_x)

#setting the parameter space
h = np.zeros([n_x, n_t])
S = np.zeros([n_x, n_t])

#solution convergence
ep_in = 0.01 #parameter for initial boundary condition solution
ep_S = 0.01 #parameter for the convergence of the final S 


#set up the solver for the boundary conditions on both h and S
R = 0.5 #surface stream
w = np.zeros([len(t)])
for j in range (len(t)):
    w[j] = 10 - 0.05*np.power(t[j], -1/5)

#INTIAL CONDITIONS 
#polynomial example to the 3rd degree
a = 0
b = -100
c = 0
d = 10
for i in range (0, n_x):
    h[i][0] = a*x[i]**3 + b*x[i]**2 + c*x[i] + d # like a linear down sloping from the initial position
    S[i][0] = (h[i][0] - h[i-1][0])/del_x #slope at (0,0) is not defined through this 
    #print(S[i][0])
S[0][0] = (h[1][0] - h[0][0])/del_x #Initialize first point

#BOUNDAY CONDITIONS
h,S = BC_Cases.setup(x, t, S, h, R, w, ep_in)

Plotter.plot_spec_t(h, x, 0)
#Plotter.plot_spec_t(h, x, int(n_t/2))
#Plotter.plot_h_t_int(h, x, t, 3)
Plotter.plot_h_t_int(h[10:], x[10:], t, 5) 
#Plotter.plot_h_t_int(S, x, t, 3)
#Plotter.plot_h_t_int(S, x, t, 5)

#may need to perform some frequency filtering on the output graphs

