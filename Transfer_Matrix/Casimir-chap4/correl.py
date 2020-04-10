#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=10
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from numpy import linalg as LA
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
def linear(x,a,b):
    return a*x+b
### Définition des matrices de transfert des différents modèles
def Trans(kbt,champ,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def mat(L):
    d=np.zeros((L+1,L+1))
    for y,i in enumerate(d):
        d[y][y] = y
    return d

def fitexp(x,a,l):
    return a*np.exp(-x/l)

### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    Z = np.sum(np.power(w,ly))

    matphi = mat(ly)

    fonction = np.empty(ly)
    for r in np.arange(ly):
        gmat = matphi*LA.matrix_power(tm,r)*matphi
        fun = 0
        for i in np.arange(ly) :
            fun += np.dot(np.dot(v[:,i],gmat),v[:,i])*np.power(w[i],ly-r)
        fonction[r] = fun/Z

    correlthermo = -1/np.log(w[-2]/w[-1])

    popt,pcov= curve_fit(fitexp,np.arange(ly),fonction,p0=[1,correlthermo])
    print(correlthermo,popt[1])
    return [correlthermo,popt[1]]


############################
# Déclaration matrices
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1/TC;

lMin = 10; lMax = 30
lSpace = np.arange(lMin,lMax)

tSpace = np.linspace(0.8,1,2)
h = 0.1
casimir = np.empty(len(lSpace))

def fitpow(x,n,a):
    return a*np.power(x,n)
def fitlin(x,n,a):
    return a*x+n

correl = np.empty((len(lSpace),2))

for i,t in enumerate(tSpace) :
    for j,L1 in enumerate(lSpace):
        f1 = intDiag(L1,1/(t*TC),h,J) 
        correl[j] = f1
    plt.plot(lSpace,correl[:,0],label=str(t))
    plt.plot(lSpace,correl[:,1],'+',label=str(t))
    print(correl[:,1])
plt.legend()
plt.show()

