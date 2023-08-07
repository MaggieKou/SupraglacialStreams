# -*- coding: utf-8 -*-
import numpy as np
import math 
"""
Created on Fri Jul 14 14:18:16 2023

function Mesher takes an arbitrary length for input vector h (at a fixed timestep) and either 
1) linearlly interpolates to fixed mesh size, nx, in the case that the list is short; or
2) condenses to fixed mesh size, nx, in the case the the list is long
This dynamically allocates but will require accompanying x values that are also interspersed similarly


function custom_mesher takes fixed x values and linear extrapolates known h and x values onto the pre-determined grid
"""
class Node:    
    def __init__(self,value):        
        self.value = value
        self.next = None

class LinkedList:
    def __init__(self,head=None):
        self.head = head    
    def append(self, new_node):
        current = self.head
        if current:
            while current.next:
                current = current.next
            current.next = new_node
        else:
            self.head = new_node
    def index(self, i):
        current = self.head
        count = 0
        value = 999999999999
        while current:
            if count+1 == i:
                value = current.value
                # print("this is the current value", current.value)
                return value 
            else:
                count+=1
                current = current.next
        return value 

    def insert(self, new_element, position):
        """Insert a new node at the given position.
        Assume the first position is "1".
        Inserting at position 3 means between
        the 2nd and 3rd elements."""
        count=1
        current = self.head
        if position == 1:
            new_element.next = self.head
            self.head = new_element
        while current:
            if count+1 == position:
                new_element.next =current.next
                current.next = new_element
                return
            else:
                count+=1
                current = current.next
            # break
        pass
    def delete(self, i):
        """Delete the indexed node."""
        count = 1
        current = self.head
        if count == i:
            self.head = current.next
        else:
            while current:
                if count == i: 
                    break
                prev = current
                current = current.next
                count+=1
            if current == None:
                return
            prev.next = current.next
            current = None
        
    def print(self):
        current = self.head
        while current:
            print(current.value)
            current = current.next
            
def mesher (h, nx):
    h_return = LinkedList(Node(h[0]))
    cur_len = len(h)
    
    for i in range(1, cur_len):
        h_return.append(Node(h[i])) #copies list over to linked list
        # print(i, h[i])
    # h_return.print()
    
    "same length"
    
    if cur_len == nx: 
        return h_return
    
    elif cur_len < nx: #shorter vector than desired mesh size
        remaining = nx - cur_len
        while remaining > 1:
            R = remaining/cur_len #number of points to be inserted between each linkedlist
            if R >=1:
                i = 1
                while i <= cur_len:
                    extrapolated = (float(h_return.index(i+1))+float(h_return.index(i)))/2 #averages adjacent points
                    # print("this is the position: ",i, "and the extrapolated figure is", extrapolated)
                    h_return.insert(Node(extrapolated), i+1) #inserts note at that position
                    i = i +2
                cur_len = 2*cur_len -1;  #update new length of linked list
                remaining = nx - cur_len 
                R = remaining/cur_len #recalculate remainder
                # print("now the length is", cur_len, "R = ", R)
                
            else:
                i = math.ceil(1/R)
                while i <= cur_len: 
                    extrapolated = (float(h_return.index(i+1))+float(h_return.index(i)))/2 #averages adjacent points
                    # print("this is the position: ",i, "and the extrapolated figure is", extrapolated)
                    h_return.insert(Node(extrapolated), i+1) #inserts note at that position
                    i = i + math.ceil(1/R)
                cur_len = 2*cur_len -1  #update new length of linked list
                remaining = nx - cur_len 
                R = remaining/cur_len #recalculate remainder
                
    else: #longer vector than desired mesh size
        remaining = cur_len -1 
        while remaining > 1:
            R = remaining/cur_len
            if R >= 1: 
                i = 1
                while i <=cur_len: 
                    h_return.delete(i+1)
                    i = i +2
                cur_len = 2*cur_len -1;  #update new length of linked list
                remaining = nx - cur_len 
                R = remaining/cur_len #recalculate remainder
                # print("now the length is", cur_len, "R = ", R)
            else: 
                i = math.ceil(1/R)
                while i <= cur_len: 
                    h_return.delete(i+1) #inserts note at that position
                    i = i + math.ceil(1/R)
                cur_len = 2*cur_len -1  #update new length of linked list
                remaining = nx - cur_len 
                R = remaining/cur_len #recalculate remainder
    
    return h_return

def custom_meshing (h, x_h, x_mesh):
    #assuming that every x_h corresponds to similarly indexed h
    h_return = np.zeros(len(x_mesh))
    j = 0
    for x in x_mesh:
        i = 0
        x_comp = x_h[i]
        while x_comp < x and i+1 < len(x_h): 
            i +=1 #loops until x_h is greater or more than the current x value 
            x_comp = x_h[i]
        if x_comp == x:
            h_return[j] = h[i]
        else: 
            if i < len(x_h): 
                h_return[j] = (h[i-1]*(x - x_h[i-1]) +h[i]*(x_h[i]-x))/(x_h[i]-x_h[i-1])
            else: 
                h_return[j] = 9999999999 #if u are getting a huge number for the h value, it was probably an indexing error
        j +=1
    return h_return
h_ll = custom_meshing([1,4,9, 16], [1,2,3,4], [1, 1.5, 2, 2.5, 3.5, 4])