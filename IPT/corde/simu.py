#!/usr/bin/python
# -*- coding: utf8
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

##############################
## Définition des fonctions ##
##############################

#Constants
from scipy import constants as CO
G = CO.value(u'standard acceleration of gravity')
M   = 1

#Définition de la force électrostatique
#Gravité
def grav():
    return [0,-M*G]

#ODE
import scipy.integrate as integrate
def equation(vec,t):
    x,v = vec[0:2],vec[2:4]
    return np.array([v,grav()]).flatten()

def contrainte(x):
    lengths = np.sqrt(np.sum(np.diff(x, axis=0)**2, axis=1)) # Length between corners
    total_length = np.sum(lengths)
    return total_lenght-MaxL

MaxL = 2
angle = 20
v0 = 1
x0 , y0 , vx0 , vy0 = 0,0,np.cos(angle)*v0,np.sin(angle)*v0
xinit = [x0,y0,vx0,vy0]

t= np.linspace(0,10,20000)
z = integrate.minimize(equation,xinit,t)

optimize.minimize(equation,init,

plt.plot(z[:,0],z[:,1])
plt.legend()
plt.axis([-1,10,-10,1])
plt.show()
