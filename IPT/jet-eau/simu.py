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

#Coefficient de frottement d'une goutte sphérique dans l'air
#F  = C rho S /2

#C  = 0.5 si Re < 3e6, C = 0.1 si Re > 3e6
#Reynolds d'une goutte dans l'air : 
#Re = VL/nu 
#   = 1*1e-2 / (1.57*1e-5)
#   = 1e3
C   = 0.5 
#rho = 1,225 kg/m^3
Rho = 1.225
#Taille moyenne d'une goute d'eau : 1e-3 m

S   = np.pi * (4e-3)**2
#
F   = C*Rho*S/2
#
#Masse d'une goutte d'eau
M   = 34e-3

#Définition de la force électrostatique
pos = [3,-4]
def elec(vec,q):
    diff = np.subtract(vec,pos)
    return q*V* diff / LA.norm(diff)**3
#Gravité
def grav():
    return [0,-M*G]

def frot(vec):
    return [-F * vec[0]**2, -F * vec[1]**2]

#ODE
import scipy.integrate as integrate
def equation(vec,t,q,k):
    x,v = vec[0:2],vec[2:4]
    return np.array([v,(grav()+elec(x,q)+frot(v)) / M]).flatten()

t = np.linspace(0,2,1e6)

x0 , y0 , vx0 , vy0 = 0,0,0,0
Q = -10
xinit = [x0,y0,vx0,vy0]

plt.plot(pos[0],pos[1],'ro')
k = 0
for V in np.linspace(0.5,1,6) :
    z = integrate.odeint(equation,xinit,t,args=(Q,k))
    plt.plot(z[:,0],z[:,1],label="V="+str(V))
plt.legend()
plt.axis([-1,10,-10,1])
plt.show()
