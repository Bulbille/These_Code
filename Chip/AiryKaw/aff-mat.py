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
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
from scipy.optimize import curve_fit
from numpy import linalg as LA

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

fichier=np.loadtxt(sys.argv[1])
T = 0.5*2/np.log(1+2**0.5)
B = 0.1
LY = len(fichier[:,0])
plt.plot(fichier[:,0],distrib(LY,1/T,B,1)

plt.plot(fichier[:,0],fichier[:,1]/np.sum(fichier[:,1]))



plt.show()

