#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=10
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import factorial
from numpy import linalg as LA
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
### Définition des matrices de transfert des différents modèles
def Trans(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def mat(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = abs(y)
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    Z = (np.sum(np.power(w,ly)))
    matphi = mat(ly)
    phi =  np.dot(np.dot(v[:,-1],matphi),v[:,-1])

    distrib = np.zeros(len(v[0]))
    for i,li in enumerate(w) :
        distrib += np.power(li,ly)*v[:,i]**2
    return distrib/Z

J = 1 ; BETA = 1 ; 
L = 40
LSpace = np.arange(L)
h = 0.01

def poisson(x,l):
    return np.power(l,x)/factorial(x)*np.exp(-l)

for L in [5,10,20,40,80] :
    distrib = intDiag(L,BETA,h,J)
    plt.plot(np.arange(L),distrib,label='L='+str(L))
#    popt,pcov = curve_fit(poisson,np.arange(L),distrib,p0=[1])
#    print(popt)
#    plt.plot(np.arange(L),poisson(np.arange(L),*popt),'+')

plt.xlabel('$p(h)$')
plt.ylabel('$h$')
plt.xlim([0,20])
#plt.yscale('log')
plt.legend()
plt.savefig('distribution-taille-finie.pdf')
plt.show()
