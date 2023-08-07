# -*- coding: utf-8 -*-

import solvers
import plotter
import numpy as np

#Initial conditions
xo = 0.
xf = 10.
del_x = 0.01
x = np.arange(xo, xf, del_x)
n_x = len(x)

to = 0.
tf = 10.
del_t = 0.01
t = np.arange(to, tf, del_t)
n_t = len(t)

v = 0.5
err = 0.7

h = np.zeros([n_x,n_t])
S = np.zeros([n_x,n_t])

#parameters for initial conditions
a = 0
b = -0.1
c = 0
d = 10
#e = -10


for i in range(n_x):
    h[i][0] = a*np.power(x[i], 3) + b*np.power(abs(x[i]), 2) + c*x[i] + d
    #h[i][0] = a*np.tanh(b*x[i] - c)+ d + e*x[i]
for i in range(1, n_x):
    S[i][0] = (h[i][0] - h[i-1][0])/del_x
S[0][0] = S[1][0] #assume relatively smooth 

h,S = solvers.init_BC(h, S, t, v)
h,S = solvers.iter_max_thresh(h, S, t, x, err)

plotter.h_spec_t(h,x,t,0)
t_steps = [0, 10, 50, 99]
plotter.h_plot_series(h,x,t, t_steps)