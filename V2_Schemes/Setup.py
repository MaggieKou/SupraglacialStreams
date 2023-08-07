# -*- coding: utf-8 -*-
"""
Defines ICs and solving 
"""
import numpy as np
import SchemeV2Solvers
import sys
sys.path.append('C:\\Users\\maggi\\Documents\\School Stuff\\ESROP\\Code\\SupraglacialStreams\\V1_Schemes')
import Plotter

"TIMESTEPS"
to = 0.1 #starting at the intial condition
tf = 100
n_t = 100 #number of timesteps 
del_t = (tf-to)/n_t
#ho = np.linspace(to, tf, n_t) #initiate initial conditions 
t = np.linspace(to, tf, n_t)

"SPATIAL STEPS"
xo = 0 #starting at the intial condition
xf = 1
n_x = 100 #number of spatial steps
del_x = (xf-xo)/n_x
#ho = np.linspace(to, tf, n_t) #initiate initial conditions 
x = np.linspace(xo, xf, n_x)

"setting the parameter space"
h = np.zeros([n_x, n_t])
S = np.zeros([n_x, n_t])

"solution convergence"
err = 0.00001 #parameter for error threshold in solver



"set up the solver for the boundary conditions on both h and S"
Q = 2 #discharge rate of the stream 
R = 0.5 #surface stream
w_t = np.zeros([len(t)])
for j in range (len(t)):
    w_t[j] = 0.000000002*t[j]**(-1000)
    #-1000*np.power(abs(S[0][0]), -1/5)/n_t #rate of increase in the width of the river, i.e. by the end it will have doubled 

"INTIAL CONDITIONS"
#polynomial example to the 3rd degree
a = 0
b = -20
c = -1
d = 100
for i in range (0, n_x):
    h[i][0] = a*x[i]**3 + b*x[i]**2 + c*x[i] + d # like a linear down sloping from the initial position
    S[i][0] = (h[i][0] - h[i-1][0])/del_x #slope at (0,0) is not defined through this 
    #print(S[i][0])
S[0][0] = (h[1][0] - h[0][0])/del_x #Initialize first point

"Setting up BCs"
S = SchemeV2Solvers.S_BC(S, w_t, del_t, Q, R)
h = SchemeV2Solvers.h_BC(del_t, h, S)

"SOLVING IN SPACETIME"
#h = SchemeV2Solvers.strang_split(h, S, n_x, n_t, del_t, del_x, x, err)
for i in range(1, n_x):
    for j in range (1, n_t):
        h[i][j], S[i][j] = SchemeV2Solvers.solve_fountain(h, S, i, j, del_t, del_x, err, 0.01)
        
for i in range(1, n_x):
    for j in range (0, n_t):
        h[i][j] = h[i-1][j] + del_x*(S[i][j]+S[i-1][j])/2
'''

'''
#h = SchemeV2Solvers.solve_spacetime(h, S, n_x, n_t, del_t, del_x, err)
#h,S = SchemeV2Solvers.smoothing(h, S)


       
#after filtering, can we feed these points back into the solver 
#incising case solver 
#Plotter.plot_spec_t(h, x, 0)
Plotter.plot_h_t_int(h[:], x[:], t, 5) 
Plotter.plot_h_t_int(S[:], x[:], t, 5) 


#may need to perform some frequency filtering on the output graphs
