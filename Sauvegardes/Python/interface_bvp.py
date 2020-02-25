#! /usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
import matplotlib.pyplot as plt 
import scipy.integrate as Int 
from scipy.integrate import solve_bvp

plt.rc("font",size=20)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)

color=['#ff0000','#ff3300','#ff6600','#ff9900','#00ff00','#009900','#006600','#003300','#00ffff','#0099ff','#0000ff','#9900ff','#cc66ff']

L   = 10

def edo(x,y,p):
    E       = p[0]
    epsilon = p[1]
    return np.vstack((y[1],-2*(E-epsilon*np.heaviside(x,0.5))*y[0]))

def bc(ya,yb,p):
    return np.array([ya[0],yb[0]])

x   = np.linspace(0,L,L*5)
y_a = np.zeros((2,x.size))
y_a[0] = 3

res = solve_bvp(edo,bc,x,y_a)

x_plot  = np.linspace(0,1,100)
y_plot  = res.sol(x_plot)[0]
plt.plot(x_plot,y_plot)
plt.show()

