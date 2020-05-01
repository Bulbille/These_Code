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

import numpy as np
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

from scipy.optimize import curve_fit
from numpy import linalg as LA
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate

### Fonction qui calcule directement l'énergie libre
def TransGauss(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*(y1-y2)**2) )
    return d

def distrib(ly,beta,h,j):
    tm = TransGauss(beta,h,j,ly)
    w,v = LA.eigh(tm)
    Z = np.sum(np.power(w,ly))
    ph = np.zeros(ly)
    for l in np.arange(ly):
        ph[l] = np.sum(np.power(w,ly)*np.power(v[l,:],2))
    return ph/Z,-1/np.log(w[-2]/w[-1])


T = 3
B = 0.01
L = 400
res = distrib(L,1/T,B,1)
print(T,B,res)
ph = res[0]
plt.plot(np.arange(L)+200,ph,'+',label='Matrice de Transfert')


glau = np.loadtxt(sys.argv[1])
kaw = np.loadtxt(sys.argv[2])

plt.plot(glau[:,0]-400,glau[:,1]/np.sum(glau[:,1]),label="Glauber")
plt.plot(kaw[:,0],kaw[:,1]/np.sum(kaw[:,1]),label="Kawasaki")

plt.xlabel('$h$')
plt.xlim([370,430])
plt.ylabel('$p(h)$')
plt.legend()
plt.savefig('comp-airy-kaw-glau.pdf')
plt.show()

