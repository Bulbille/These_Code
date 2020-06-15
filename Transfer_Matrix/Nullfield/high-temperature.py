#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=16
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import brentq
from numpy import linalg as LA
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
### Définition des matrices de transfert des différents modèles
def Trans(kbt,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*inter*abs(y1-y2))
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,j):
    tm = Trans(beta,j,ly)
    w,v = LA.eigh(tm) ;
    return -np.log(w[-1])/beta

#def equaDiag(ly,beta,j):
#    betexp = lambda beta : np.exp(-beta)
#    theta = lambda eig,h,r : 1/np.tan(h*eig) - r*np.sin(eig)/(1-r*np.cos(eig))
#    ret:

############################
# Déclaration matrices
J = 1
LY = 100

def dean(j,ly,beta):
    return -np.log(ly+1-beta*j*(ly**2+2*ly)/3)/beta

bs = np.logspace(-4,-1,60)
l0 = np.empty(len(bs))
l1 = np.empty(len(bs))

for i,BETA in enumerate(bs) :
    res = intDiag(LY,BETA,J)
    l0[i] = res

plt.plot(bs,l0,label='$f(\\beta,L=100)$')
plt.plot(bs,dean(J,LY,bs),'+',label='$- \\frac{1}{\\beta} \log(L+1-\\beta J (L^2+2L)/3)$')
#print(dean(J,LY,bs))


plt.legend(loc='lower right')
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('$\\beta$')
plt.ylabel('$f(\\beta,L=100)$')

plt.tight_layout()
plt.savefig('high_temperature.pdf')
plt.show()
