#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 14:35:32 2023

@author: mk2217
"""

import ODE_solve
import plotter
import numpy as np

#Initial conditions
xo = 0.
<<<<<<< HEAD
xf = 5.
del_x = 0.01
=======
xf = 10.
del_x = 1
>>>>>>> 160244ca137cc412391c4e8781a45791bb2d755f
x = np.arange(xo, xf, del_x)
n_x = len(x)

to = 0.
tf = 10.
<<<<<<< HEAD
del_t = 0.1
t = np.arange(to, tf, del_t)
n_t = len(t)

v = 0.7

H_in = np.zeros([n_x])

#polynomial function
a = -2
b = 2
c = 4
d = 10
e = -1

'''
#tanh 
a = -3
b = 1
c = 2
d = 9
'''

for i in range(n_x):
    #H_in[i] = a*np.power(x[i], 3) + b*np.power(abs(x[i]), 2) + c*x[i] + d
    H_in[i] = a*np.tanh(b*x[i] - c)+ d + e*x[i]

h, s = ODE_solve.timestepper(H_in, x, t, v)
plotter.h_spec_t(h,x,t,0)
t_steps = [0, 25, 50, 99]
=======
del_t = 1
t = np.arange(to, tf, del_t)
n_t = len(t)

v = 10

H_in = np.zeros([n_x])
#tanh 

#polynomial function
a = 0
b = 0
c = -0.1
d = 1000

for i in range(n_x):
    H_in[i] = a*np.power(x[i], 3) + b*np.power(abs(x[i]), 2) + c*x[i] + d

h, s = ODE_solve.timestepper(H_in, x, t, v)
plotter.h_spec_t(h,x,t,0)
t_steps = [0, 2]
>>>>>>> 160244ca137cc412391c4e8781a45791bb2d755f
plotter.h_plot_series(h,x,t, t_steps)