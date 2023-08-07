#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 14:42:36 2023

@author: mk2217
"""
from matplotlib import pyplot as plt
import numpy as np

def h_space(h, x, t_j):
    n_x = len(h)
    h_ints = np.zeros([len(x)])
    
    for i in range(0, n_x):
        h_ints[i] = h[i][t_j]
    return h_ints

def h_spec_t (h,x,t, t_j):
    h_ints = h_space(h,x,t_j)
    fig, ax = plt.subplots()
    h_plt = ax.plot(x, h_ints)
    plt.xlabel('x [nondim]')
    plt.ylabel('h [nondim]')
    plt.title('Solution to h at t ='+ str(t[t_j]))

def h_plot_series(h, x, t, t_steps):
    times = np.zeros([len(t_steps)])
    h_ser = np.zeros([len(t_steps),len(x)])
    for j in range(0, len(t_steps)):
        times[j] = str(np.round(t[t_steps[j]], 1)) #labels times
        h_ser[j] = h_space(h, x, t_steps[j])
    fig, ax = plt.subplots()
    for j in range(0, len(t_steps)):
        h_plot = ax.plot(x, h_ser[j])
    
    ax.legend(times)
    plt.xlabel('x [nondim]')
    plt.ylabel('h [nondim]')
    #ax.set_yscale("log", base = 10)
    plt.title('Solution to h at various times')