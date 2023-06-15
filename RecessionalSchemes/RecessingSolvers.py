# -*- coding: utf-8 -*-
"""
SOLVERS and HELPING FUNCTIONS FOR THE RECEEDING PROFILE CASE,

given the initial condition H'(0) - this will give a stead state spatial profile
our parameterization of the profile in 3D can be solved for with the additional v dimension
"""

import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
"helper functions"

#these are both 1D vectors
def small_slope(H_p):
    H_pp = -5/2*np.sign(H_p)*np.power(abs(H_p), 11/5)*(np.sign(H_p)*np.power(abs(H_p), 1/5) + 1)
    return H_pp

def large_slope(H_p):
    H_pp = -5/2*np.sign(H_p)*np.power(abs(H_p), 11/5)
    return H_pp

def general_slope(v, x, S, H):
    del_x = (x[len(x)-1]-x[0])/len(x)
    for i in range(len(v)):
        for j in range (1, len(x)):
            H_pp = -5/2*np.sign(S[i][j-1])*np.power(abs(S[i][j-1]), 11/5)*(np.sign(S[i][j-1])*np.power(abs(S[i][j-1]), 1/5)/(np.sqrt(1+np.power(v[i], 10)*np.power(abs(S[i][j-1]), 2))) + 1)
            print(H_pp)
            S[i][j] = S[i][j-1] + del_x*H_pp
            H[i][j] = H[i][j-1] + del_x*(S[i][j-1] + S[i][j])/2
            
    return S, H
#parameters for the recessional velocity, if this is a constant. 
n_v = 5 #steps for v
v_lower = 0.45
v_upper = 5
v = np.linspace(v_lower, v_upper, n_v)
H_p = -0.1 #initial slope
H_in = 1000 #initial height of the stream

n_x = 10000
xo = 0
xf = 80
x = np.linspace(xo, xf, n_x)
H = np.zeros([n_v, n_x])

S_out = np.zeros([n_v, n_x]) # this reaches a singularity lol
for i in range (n_v):
    S_out[i][0] = H_p #IC set
    H[i][0] = H_in #IC for position 

S_out, H = general_slope(v, x, S_out, H)

def plot_H(H,x, v):
    fig,ax = plt.subplots()
    time_label = len(v)*["a"]
    for j in range(len(v)):
        h_plt = ax.plot(x, H[j]) #all the times plotted
        time_label[j] = "v = " + str(np.round(v[j], 2))
        #print(h_ints)
    ax.legend(time_label) #there may be an issue with global permissions for this to show up 
    ax.set_yscale("log", base = 10)
    #ax.set_xscale("log")
    plt.xlabel("x [dimensionless]")
    plt.ylabel("H [dimensionless]")
    plt.title("Solving Scheme on H for different receding velocities")
    plt.show()
        
""" COLOUR MAP
fig, ax = plt.subplots(1, 1)
pcm = ax.pcolor(v[1:], x, H[1:][:],
                   norm=colors.Normalize(vmin=H[1:][:].min(), vmax=H[1:][:].max()),
                   cmap='viridis', shading='nearest')
fig.colorbar(pcm, ax = ax, extend='max')


#fig.colorbar(im, ax=ax, label='Interactive colorbar')
"""
plot_H(H, x, v)