# -*- coding: utf-8 -*-
"""
Created on Tue Jul 08 20:43:22 2014

@author: gasmgpu1
"""

import numpy as np
# NORM OF A 3D VECTOR
def nrm(v):
    return (v[0]**2 + v[1]**2 + v[2]**2)**0.5    

# CROSS PRODUCT OF 3D VECTORS
def crs(v1, v2):
    x = ((v1[1] * v2[2]) - (v1[2] * v2[1]))
    y = ((v1[2] * v2[0]) - (v1[0] * v2[2]))
    z = ((v1[0] * v2[1]) - (v1[1] * v2[0]))
    v = np.array( [x,y,z] )
    return v

# DOT PRODUCT OF 3D VECTORS    
def dot(v1, v2):
    s = (v1[0] * v2[0]) + (v1[1] * v2[1])  + (v1[2] * v2[2]) 
    return s

# 3x3 IDENTITY MATRIX    
def eye():    
    eye = np.array( [[1.0, 0.0, 0.0],[0.0, 1.0, 0.0],[0.0, 0.0, 1.0]] )
    return eye    
def dia(n):    
    dia = np.array( [[n, 0.0, 0.0],[0.0, n, 0.0],[0.0, 0.0, n]] )
    return dia 

def do3(m1, m2, m3):
    return np.dot(m1,np.dot(m2, m3))
# ANTICLOCKWISE NORMAL OF SEGMENTS   
def normal(p1, p2):
    dp =  p2 - p1
    sv = dp / nrm( dp )     # versor tangente 
    nvt = np.array( [ 0.0, sv[2], -sv[1] ])
    return (nvt/nrm(nvt))
    
def tangen(p1, p2):
    tv = (p2 - p1)
    return (tv/nrm(tv))

def versor(v):    
    return  (v / ( (v[0]**2 + v[1]**2 + v[2]**2)**0.5 )) 

def unitx():
    return np.array( [1.0, 0.0, 0.0] )

def unity():
    return np.array( [0.0, 1.0, 0.0] )

def unitz():
    return np.array( [0.0, 0.0, 1.0] )