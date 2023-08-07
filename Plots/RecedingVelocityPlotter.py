# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt

"""
Created on Thu Jul 13 09:26:48 2023

@author: maggi
"""
x = np.linspace(0.95, 1, 1000)


def exp_ret (x, a, b, c, d, hwhm):
    ret = a*np.exp(b*x) + c*np.exp(d*x)
    for i in range(0, len(ret)):
        if x[i] < 0.5 + hwhm:
            ret[i] = 0
    return ret

v_vec = np.zeros([9, len(x)])
v_vec[0] = exp_ret(x, 0.4944, -8.9275, 0, -7.2569, 0.1101825)
v_vec[1] = exp_ret(x, 0.3462, -8.1684, 0.000, 0.9112, 0.1101875)
v_vec[2] = exp_ret(x, 0.5517, -8.6735, 0.000, -0.2899, 0.110194)
v_vec[3] = exp_ret(x, 0.6119, -8.8050, 0, 0.000, 0.110227)
v_vec[4] = exp_ret(x, 0.5572, -8.8146, 0, 0, 0.110282)
v_vec[5] = exp_ret(x, 0.3790, -8.3865, 0.0002, -2.1427, 0.1103295)
v_vec[6] = exp_ret(x, 0.2494, -8.3002, 0.032, -6.6944, 0.110393)
v_vec[7] = exp_ret(x, 0.3806, -8.2165, -1.2540, -11.8150, 0.110725)
v_vec[8] = exp_ret(x, 0.335, -8.2675, -13.4647, -15.7317, 0.11128)

fig, ax = plt.subplots()

for i in range(0, 9):
    h_plot = ax.plot(x, v_vec[i])
   
ax.legend(["hx = 10^-3", "hx = 7*10^-4", "hx = 5*10^-4", "hx = 2*10^-4", "hx = 10^-4","hx = 7*10^-5", "hx = 5*10^-5", "hx = 2*10^-5", "hx = 10^-5"])
plt.xlabel('x [nondim]')
plt.ylabel('v_rec [nondim]')
    #ax.set_yscale("log", base = 10)
plt.title('Receding velocity as a function of position at different initial slopes')
