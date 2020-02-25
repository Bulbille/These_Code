#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=9
labelsize=16
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
from numpy import linalg as LA
#np.set_printoptions(precision=2)
import math as mt
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

TC = 2/np.log(1.+np.power(2,0.5));BETA = 1/TC
LY = 400; J = 1;  

colors=['red','blue','green','rosybrown','sandybrown','indigo','darkmagenta','gold','olivedrab'] ; nb_colors = len(colors)

### Définition des matrices de transfert des différents modèles
def lnfact(n):
#    return mt.log(mt.factorial(n))
    if(n<25):
        return mt.log(mt.factorial(n))
    else:
        return n*mt.log(n)-n
def TransPOP(kbt,inter,L,mu):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*(-mu*(y1+y2)/2+inter*abs(y1-y2))-(lnfact(y1)+lnfact(y2))/2 )
    return d

### Définition des matrices de magnetisation des différents modèles
def Mat(L):
    d=np.zeros((L+1,L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(y)
    return d

### Fonction qui calcule directement <H>, <H^2> et E

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


varplot = plt.subplot(211)
hplot = varplot.twinx()
eplot = plt.subplot(212)

nbmu = 50
mumax = 20
mubar = np.linspace(0,mumax,nbmu)
res = np.empty((nbmu,4))
for nb,MU in enumerate(mubar):
    diag = intDiag(LY,BETA,J,MU)
    res[nb,0] = MU
    res[nb,1] = diag[0]
    res[nb,2] = diag[1]
    res[nb,3] = diag[3]
    print(MU,diag[0])

hplot.plot(res[:,0],res[:,1],color='red')
varplot.plot(res[:,0],res[:,2],color='blue',label="$<H^2>$")
eplot.plot(res[:,0],res[:,3])


varplot.set_xlabel('$f$')
varplot.set_ylabel('$\sigma = (<H^2>-<H>^2)^{0.5}$',fontsize=labelsize)
varplot.legend(loc='upper left')

varplot.set_ylabel('$<H>$',fontsize=labelsize)
hplot.legend(loc='upper right')

varplot.set_xlabel('$f$')
varplot.set_ylabel('$E = - \\frac{\partial \ln(Z)}{\partial \Beta}$',fontsize=labelsize)
varplot.legend(loc='upper left')


np.savetxt('calculmat',res)
plt.show()
plt.close()


