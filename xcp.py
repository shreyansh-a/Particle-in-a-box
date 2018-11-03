# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 10:54:36 2018

@author: nitin
"""
import numpy as np
from fractions import Fraction
from scipy.linalg import eigh
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib as mpl

def get_PZxcpot(d):
    d=d.flatten()
    n2=int(len(d))
    Vx=np.zeros((n2,n2))
    rs=np.zeros(n2)
    ec=Vx
    rsecdr=Vx
    
    for i in range(n2) :
        Vx[i][i] = ((3/np.pi)*d[i])**Fraction('1/3')
        "PERDEW-ZUNGER CORRELATION"
        rs[i] = (3/(4*np.pi*d[i]))**Fraction('1/3')
        if rs[i] < 1 :
            ec[i][i] = -0.0480+0.031*np.log(rs[i])-0.0116*rs[i]+0.0020*rs[i]*np.log(rs[i])
            rsecdr[i][i] = -rs[i]*((0.031/rs[i])-0.0116+0.0020*(1+np.log(rs[i])))/3
            
        else :
            ec[i][i] = -0.1423/(1+1.0529*np.sqrt(rs[i])+0.3334*rs[i])
            rsecdr[i][i] = (-0.1423*(1.0529*np.sqrt(rs[i])+2*0.3334*rs[i]))/(6*(1+1.0529*np.sqrt(rs[i])+0.3334*rs[i])**2)
    vc= ec+rsecdr
    V=Vx+vc
    return V
  

def get_den(x1,x2,y1,y2,z1,z2,steps):
    x=np.linspace(x1,x2,steps)
    y=np.linspace(y1,y2,steps)
    z=np.linspace(z1,z2,steps)
    g=np.zeros(steps**3)
    m=0
    for i in range(steps):
        for j in range(steps):
            for k in range(steps):
                    g[m]=eval("np.exp(x[k]**2+y[j]**2+z[i]**2)")
                    m=m+1
    return g
    
n=get_den(0,.5,0,.6,0,.7,10)
v=get_PZxcpot(n)