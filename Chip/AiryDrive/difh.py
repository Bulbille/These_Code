#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=14
labelsize=16
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
from numpy import linalg as LA
import math as mt
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

BETA = 1
LY = 100; J = 1;

def Trans(kbt,inter,L,mu):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*(champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*abs(y1-y2)))
    return d

### Définition des matrices de magnetisation des différents modèles
def Mat(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = abs(L/2-y)
    return d

### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,j,mu):
    matphi = Mat(LY)
    matphi2 = np.square(matphi)
    matphi3 = np.power(matphi,3)

    tm = TransPOP(beta,j,ly,mu)
    w,v = LA.eigh(tm) 

    phi =  np.dot(np.dot(v[:,-1],matphi),v[:,-1])
    phi2 =  np.dot(np.dot(v[:,-1],matphi2),v[:,-1])
    sigma = (phi2-phi**2)**0.5
    gamma = 0 

    #Dérivée à l'ordre un avec 5 points E = -d(lnZ)/dB
    # Finite difference calculator http://web.media.mit.edu/~crtaylor/calculator.html
    db = 0.001
    i=2 
    f = np.empty(2*i+1)
    for n in np.arange(2*i+1):
        tm = TransPOP(beta+(n-i)*db,j,ly,mu)
        w,v = LA.eigh(tm) 
        f[n] = mt.log(w[-1])
    E =  (1*f[i-2]-8*f[i-1]+0*f[i+0]+8*f[i+1]-1*f[i+2])/(12*1.0*db**1)
    E = mu*phi - E 

    return [phi,sigma,gamma,E]


###############################
# Mise en mémoire des données #
###############################

energie = plt.subplot(111)
derivene = energie.twinx()

data = np.loadtxt(sys.argv[1])

energie.plot(data[:,0],data[:,1],color="blue")
derivene.plot(data[:,0],np.gradient(data[:,1],data[:,0]),color="green")

energie.set_xlabel('$f$')
derivene.set_xlabel('$\\frac{d E}{df}$')
energie.set_ylabel('E')

plt.savefig('shear-ene-sos.pdf')
plt.show()
plt.close()
