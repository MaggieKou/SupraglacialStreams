# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 14:01:44 2023

@author: maggi
"""

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm

#produces the phase plot of H' for varying values of v 

H_p_o = -2
H_p_f = 2
n_H = 1000
H_p = np.linspace(H_p_o, H_p_f, n_H)

n_v = 1000
v_o = 0.7
v_f = 1.2
v = np.linspace(v_o, v_f, n_v)
H_pp_mult = np.zeros([n_v, n_H])

for i in range (0, n_v):
    for j in range (0, n_H):
        H_pp_mult[i][j] = -5/2*np.sign(H_p[j])*np.power(abs(H_p[j]), 11/5)*(np.sign(H_p[j])*np.power(abs(H_p[j]), 1/5)/(np.sqrt(1+np.power(v[i], 10)*np.power(abs(H_p[j]), 2))) + 1)


#colourmap plotter 

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(v, H_p, H_pp_mult, vmin=H_pp_mult[:][:].min(), cmap=cm.plasma)

ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[], xlabel = "v", ylabel = "H'", zlabel = "H''")

plt.show()




#fig.colorbar(im, ax=ax, label='Interactive colorbar'