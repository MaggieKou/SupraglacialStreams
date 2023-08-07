# -*- coding: utf-8 -*-

import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import cm

#produces the phase plot of H' for varying values of v 

S_o = -5
S_f = 5
n_S = 1000
S = np.linspace(S_o, S_f, n_S)

n_ht = 1000
ht_o = -2
ht_f = 0
ht = np.linspace(ht_o, ht_f, n_ht)

St = np.zeros([n_S, n_ht])
St2 = np.zeros([n_S, n_ht])

for i in range (0, n_S):
    for j in range (0, n_ht):
        St[i][j] = -5/2*np.power(abs(S[i]), 6/5)*((np.power(abs(S[i]), 6/5))/(np.sqrt(1+np.power(abs(S[i]),2))) + ht[j])
        St2[j][i] = -5/2*np.power(abs(S[i]), 6/5)*((np.power(abs(S[i]), 6/5))/(np.sqrt(1+np.power(abs(S[i]),2))) + ht[j])
#colourmap plotter 

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.plot_surface(ht, S, St2, vmin=St.min(), cmap=cm.plasma)

ax.set(xticklabels=[],
       yticklabels=[],
       zticklabels=[], xlabel = "ht", ylabel = "S", zlabel = "St")


plt.show()




#fig.colorbar(im, ax=ax, label='Interactive colorbar'