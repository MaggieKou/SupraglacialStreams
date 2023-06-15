# -*- coding: utf-8 -*-
"""
Created on Tue May 30 15:29:25 2023

@author: maggi
"""

from matplotlib import pyplot as plt 
import numpy as np

def plot_spec_t(h,x,t_j):
    n_x = len(h)
    h_ints = np.zeros([len(x)])

    for i in range (0, n_x):
        h_ints[i] = h[i][t_j]
    fig,ax = plt.subplots()
    h_plt = ax.plot(x, h_ints)
    plt.xlabel("x [dimensionless]")
    plt.ylabel("h [dimensionless]")
    plt.title("Solving Scheme on h at discrete time")
    plt.show()
   
    
def plot_h_t_int(h, x, t, t_intervals):
    n_x = len(h)
    n_t = len(h[0])
    h_ints = np.zeros([t_intervals, len(x)])
    times = np.zeros([t_intervals])
    time_label = t_intervals*["a"]
    for i in range (0, n_x): 
        for j in range (0, len(h_ints)):
            h_ints[j][i] = h[i][int(n_t/(t_intervals - j))-1] #samples the discrete times
            time_label[j] = "t = " + str(round(t[int(n_t/(t_intervals - j))-1],3))
    #print(h_ints)
    fig,ax = plt.subplots()
    for plt_time in range (0, len(h_ints)):
        h_plt = ax.plot(x, h_ints[plt_time]) #all the times plotted
  
    ax.legend(time_label) 

    plt.xlabel("x [dimensionless]")
    plt.ylabel("h [dimensionless]")
    plt.title("Solving Scheme on h for different times")
    plt.show()
   
    
    

