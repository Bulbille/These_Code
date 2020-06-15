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
    return w[-1],w[-2]

#def equaDiag(ly,beta,j):
#    betexp = lambda beta : np.exp(-beta)
#    theta = lambda eig,h,r : 1/np.tan(h*eig) - r*np.sin(eig)/(1-r*np.cos(eig))
#    ret:

############################
# Déclaration matrices
J = 1
BETA = 1

def dean(j,ly,n):
    return np.sinh(j)/(np.cosh(j)-np.cos(n*np.pi/ly))

ls = np.arange(100)+3
l0 = np.empty(len(ls))
l1 = np.empty(len(ls))

for i,LY in enumerate(ls) :
    res = intDiag(LY,BETA,J)
    l0[i] = res[0]
    l1[i] = res[1]

plt.plot(ls,l0,label='$\lambda_0$')
plt.plot(ls,dean(BETA*J,ls,1),label='$\\frac{\sinh(\\beta J) }{\cosh(\\beta J)-\cos(\pi/L)}$')
plt.plot(ls,l1,label='$\lambda_1$')
plt.plot(ls,dean(BETA*J,ls,2),label='$\\frac{\sinh(\\beta J) }{\cosh(\\beta J)-\cos(2 \pi/L)}$')
#plt.plot(ls[::5],dean(BETA*J,ls[::5],0),'+-',color="black",label='$\lambda_0(L\\to \infty)$')


plt.legend()
plt.xlabel('$L$')
plt.ylabel('$\lambda$')
plt.tight_layout()

plt.savefig('null_deanJ.pdf')
plt.show()
