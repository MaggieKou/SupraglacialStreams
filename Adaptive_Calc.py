# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 16:15:50 2023

adaptive size x mesh, all calculated for a specific point
"""
import sys
import Mesher
import numpy as np
from matplotlib import pyplot as plt
sys.path.append('C:\\Users\\maggi\\Documents\\School Stuff\\ESROP\\Code\\SupraglacialStreams\\V3_Schemes')
import plotter
import time


'''rewrite thresh_calc to include upper threshold PLS'''
def thresh_calc(dt, Q, thresh_percent, hx, h_p, hx_p, h_o, err):
    #amount to decrease the precision by to get a better mesh, want to increase the spatial mesh
    
    #given fixed dt, inputted values of Q and h_x, choose mesh size based on percent of threshold
    # SMALL SLOPE: dx = thresh_percent*dt*6/5*np.power(Q, 3/5)*np.power(abs(hx), 1/5)
    dx = dt*6/5*np.power(Q, 3/5)*np.power(abs(hx), 1/5)/np.power(abs(1 + np.power(hx,2)), 8/5)
    
    #recalculate temp values for iteration
    hx_temp = hx_iter(h_o, h_p, hx_p, dx, dt, Q, err)
    v_temp = dx/dt #value given present 
    v_prop_temp = 6/5*np.power(Q, 3/5)*np.power(abs(hx_temp), 1/5)/np.power(abs(1 + np.power(hx_temp,2)), 8/5)
    # small slope 6/5*np.power(Q, 3/5)*np.power(hx_temp, 1/5) #not thresholded
    #print("this is the propagation velocity", v_prop_temp, "this is the dx/dt", dx/dt, "the difference is", dx/dt- v_prop_temp)
    while v_prop_temp > v_temp:
        dx = dt*thresh_percent*6/5*np.power(Q, 3/5)*np.power(hx_temp, 1/5)/np.power(abs(1 + np.power(hx_temp,2)), 8/5)
        #print("dx", dx)
        #recalculate temp values for iteration
        hx_temp = hx_iter(h_o, h_p, hx_p, dx, dt, Q, err)
        
        # if hx_temp > 1: #if slope surpasses threshold slope of 1, we jump to the end of the mesh 
        #     dx = 100*dx
        #     hx_temp = hx_iter(h_o, h_p, hx_p, dx, dt, Q, err)
            
        v_temp = dx/dt
        v_prop_temp = 6/5*np.power(Q, 3/5)*np.power(abs(hx_temp), 1/5)/np.power(abs(1 + np.power(hx_temp,2)), 8/5)
        thresh_percent = thresh_percent + 0.01 #add one percent to the increased mesh value corresponding to 10 m
       
        #small slope: v_prop_temp = 6/5*np.power(Q, 3/5)*np.power(hx_temp, 1/5)
       # print("this is the propagation velocity", v_prop_temp, "this is the dx/dt", dx/dt, "the difference is", dx/dt- v_prop_temp)
       
       #threshold for the v_prop term 
    return dx

def hx_iter(h_o, h_p, hx_p, dx, dt, Q, err):
    #calculates values for hx from h_p, converge y integrating to find ht
    #h_p is guess for previous h? 
    dhdt = (h_o-h_p)/dt #guess same as previous point
    hx = hx_p + dt*5/2*(np.power(Q, -2/5)*np.power(abs(hx_p), 6/5))*(dhdt - np.power(Q, 3/5)*np.power(abs(hx_p), 6/5))
    
    #new point 
    h_temp = h_o + dx*hx
    dhdt = (h_temp-h_p)/dt
    hx_iter = hx_p + dt*5/2*(np.power(Q, -2/5)*np.power(abs(hx_p), 6/5))*(dhdt - np.power(Q, 3/5)*np.power(abs(hx_p), 6/5))
    
    while(abs(hx_iter - hx)> err):
        hx = hx_iter
        h_temp = h_o + dx*hx
        dhdt = (h_temp-h_p)/dt
        hx_iter = hx_p + dt*5/2*(np.power(Q, -2/5)*np.power(hx_p, 6/5))*(dhdt - np.power(Q, 3/5)*np.power(hx_p, 6/5))
    return hx_iter

#solves for a given timestep with initial data from the previous timestep 
def t_solve(h_o, x_p, h_p, hx_p, dt, Q, err):
    #h_o is BC for the given time step, h_p is the previous time step (vector), x_p is the previous timestep's grid 
    nx = len(x_p)
    dx = (x_p[nx-1] - x_p[0])/nx
    h_out = Mesher.LinkedList(Mesher.Node(h_o)) #used linkedlist for the first point 
    hx_out = Mesher.LinkedList(Mesher.Node(0))
    x_out = Mesher.LinkedList(Mesher.Node(x_p[0])) #linkedlist for the first point
    
    i = 1
    x_out.append(Mesher.Node(float(x_out.index(i) + dx)))
    x_out_cur = x_out.index(i+1) #current point, with dx guessed
    xf = x_p[nx-1] #final point 
    while x_out_cur < xf: #recall that linked list is indexed higher than list
    #linearly extrapolate the corresponding h value here
        for j in range (1, nx): 
            x_p_comp = x_p[j-1] #value to compare in previous mesh
            x_p_comp2 = x_p[j] #value in next point
            x_comp = x_out_cur
            if x_p_comp > x_comp and x_out.index(i) < x_p_comp2: #if between two previous mesh points, then linearly extrapolate 
                h_p_temp = (h_p[j-1]*(x_out_cur-x_p[j-1]) + h_p[j]*(x_p[j] - x_out_cur))/(x_p[j]-x_p[j-1])
                hx_p_temp = (hx_p[j-1]*(x_out_cur - x_p[j-1]) + hx_p[j]*(x_p[j] - x_out_cur))/(x_p[j]-x_p[j-1])
                break 
            elif x_p_comp < x_comp and x_out.index(i) < x_p_comp2: #between current x and nsh pt from previous timestep 
                h_p_temp = (h_out.index(i)*(x_out_cur - x_comp) + h_p[j]*(x_p[j] - x_out_cur))/(x_p[j]-x_comp)
                hx_p_temp = (h_out.index(i)*(x_out_cur - x_comp) + hx_p[j]*(x_p[j] - x_out_cur))/(x_p[j]-x_comp)
                break 
            elif x_p_comp == x_out.index(i):
                h_p_temp = h_p(j-1)
                hx_p_temp = hx_p(j-1) #meshing points are the same
                break 
  
    #find temporary value 
        hx_temp = hx_iter(float(h_out.index(i)),h_p_temp, hx_p_temp, dx, dt, Q, err)
        
    #calculates threshold meshing value, and recalculates the hx value 
        dx = thresh_calc(dt, Q, 1.0, hx_temp, h_p_temp, hx_p_temp, float(h_out.index(i)), err) 
        hx_temp = hx_iter(float(h_out.index(i)),h_p_temp, hx_p_temp, dx, dt, Q, err)
        
    #appends the new, iterated value for h
        hx_out.append(Mesher.Node(float(hx_temp)))
        h_out.append(Mesher.Node(float(h_out.index(i))+ dx*hx_temp)) 
        x_out.append(Mesher.Node(float(x_out.index(i) + dx)))
        i+=1 #continues
        x_out_cur = float(x_out.index(i+1) + dx)
        #print("previous x", x_out.index(i), "updated point - x:", x_out_cur, "dx:", dx, "h_cur:", h_out.index(i), "hx:", hx_out.index(i))

        #benchmark speed here 
    
    #return three linked lists for x, h, h_x

    return x_out, h_out, hx_out #returns slope and h at each point in original mesh

def BC_solver(h_x, h, t):
    #given vector for length of t, initial height, h, and intitial slope h_x, return solved h at the next timestep
    #recall that H^2_eta = 5/2*Q^(-2/5)*H_eta*(1-H_eta^(6/5)) - 1/H_eta
    
    pass
def linked_to_list(LL):
    cur = LL.head
    counter = 0
    while cur.next != None:
        counter +=1
        cur = cur.next
    list_out = np.zeros(counter)
    
    for i in range(0, counter):
        list_out[i] = LL.index(i+1)
    return list_out

def plot_steps(h, x, t_j, t):
    #h and x are lists of LL, t is a list of indexed times to print at
    fig, ax = plt.subplots()
    #plt.hold(True)
    times = np.zeros(len(t_j))
    count = 0
    for i in t_j: 
        h_temp = linked_to_list(h[i])
        x_temp = linked_to_list(x[i])
        times[count] = str(np.round(t[i], 1)) #labels times
        if len(x_temp) == len(h_temp):
            h_plt = ax.plot(x_temp, h_temp)
        elif len(x_temp) == len(h_temp) +1:
            h_plt = ax.plot(x_temp[0:-1], h_temp)
        count += 1
    plt.title('h for a vairety of times')
    ax.legend(times)
    plt.xlabel('x [nondim]')
    plt.ylabel('h [nondim]')
    plt.show()
    pass

def curve_fitter():
    #function to smooth out solutions from linear extrapolation from before
    
    #for each couple of x values, fit to a specific curve? 
    pass


def ht_rec_calc(hx, Q):
    return 1e6*np.power(hx, 22/5)*np.power(Q, 3/5)/(2/5*np.power(Q, 2/5) +np.power(hx,16/5))
"""
Initial Conditions and Solving Parameters
"""

#MESHING
nx = 100
xo = 0
xf = 1
x = np.linspace(xo, xf, nx)
dx = (xf-xo)/nx #temporary, evenly spaced

nt = 50
to = 0
tf = 50
t = np.linspace(to, tf, nt)
dt = (tf-to)/nt

#tanh settings
a = 1e-3
b = 8
c = 0.5
d = 1e-3
e = 1e-3

#polynomial settings
# a = 1e-6
# b = 1e-4
# c = 1e-3
# d = 0

Q = 5e-3
err = 0.001

#linkedlists for initial conditions
h_init  = Mesher.LinkedList(Mesher.Node(a*np.tanh(b*(x[0] - c)) + d*x[0] + e)) #used linkedlist for the first point 
hx_init = Mesher.LinkedList(Mesher.Node(a*b*np.power(np.cosh(b*(x[0] - c)),-2) + d)) #used linkedlist for the first point 
x_init = Mesher.LinkedList(Mesher.Node(x[0]))
for i in range(1, nx):
    h_init.append(Mesher.Node(a*np.tanh(b*(x[i] - c)) + d*x[i] + e)) #append initial condition points to linkedlist 
    hx_init.append(Mesher.Node((h_init.index(i+1)-h_init.index(i))/dx))
    x_init.append(Mesher.Node(x[i]))



x_mesh = nt*[Mesher.LinkedList(Mesher.Node(0))]
h = nt*[Mesher.LinkedList(Mesher.Node(0))]
h_x = nt*[Mesher.LinkedList(Mesher.Node(0))]

"list not accepting LL datatype as argument, wants float"
h[0] = h_init
x_mesh[0] = x_init
h_x[0] = hx_init
h_xt = 5e-14 + 2.45e-11*t

h_temp_next = Mesher.LinkedList(Mesher.Node(h_init.index(1) + dt*(ht_rec_calc(hx_init.index(1), Q))))
#h_temp_next = Mesher.LinkedList(Mesher.Node(h_init.index(1) + dt*(2/5*h_xt[0]*np.power(h_init.index(1),(-6/5))*np.power(Q, (2/5)) + np.power(Q, (3/5))*np.power(h_init.index(1),(6/5)))))
h[1] = h_temp_next

temp_time = time.time()

for j in range (1, nt):
    temp_time = time.time()
    x_mesh[j], h[j], h_x[j] = t_solve(h[j].index(1), linked_to_list(x_mesh[j-1]), linked_to_list(h[j-1]), linked_to_list(h_x[j-1]), dt, Q, err)
    if j < nt- 1: 
        #print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAa")
        #h_temp_next = Mesher.LinkedList(Mesher.Node(h[j].index(1) + dt*(2/5*h_xt[j]*np.power(h[j].index(1),(-6/5))*np.power(Q, (2/5))+np.power(Q, (3/5))*np.power(h[j].index(1),(6/5)))))
        h_temp_next = Mesher.LinkedList(Mesher.Node(h[j].index(1) + dt*((ht_rec_calc(h_x[j].index(2), Q)))))
        h[j+1] = h_temp_next
    print("timestep", j, "BC: h=",h[j].index(1), "S = ", h_x[j-1].index(2), "time passed in seconds:", time.time() - temp_time)
    temp_time = time.time()
        
        #print(h[j].index(1) + dt*(np.power(Q, (3/5))*np.power(h[j].index(1),(6/5))))
# x_rec, v_rec = receding_rate(h_x, x, t)
# fig, ax = plt.subplots()
# h_plt = ax.plot(x, h_x[0][:], x, h_x[50][:],x, h_x[100][:],x, h_x[150][:],x, h_x[169][:], marker = ".")
# plt.title('Solution to h_x with adaptive solver')
# ax.legend(["t = 85", "t = 87.5", "t = 89.5", "t = 90", "t=95"])
# plt.xlabel('x [nondim]')
# plt.ylabel('hx [nondim]')
x_print = np.zeros([nt, nx])
h_print = np.zeros([nt, nx])
hx_print = np.zeros([nt, nx])


# fig, ax = plt.subplots()
# h_plt = ax.plot(t, x_rec, marker = ".")
# plt.title('position as a function of time for max sloep thresh')
# plt.xlabel('t [nondim]')
# plt.ylabel('x [nondim]')


# fig, ax = plt.subplots()
# h_plt = ax.plot(x_rec, v_rec, marker = ".")
# plt.title('Receding veolcity as a function of position')
# plt.xlabel('x [nondim]')
# plt.ylabel('v [nondim]')
