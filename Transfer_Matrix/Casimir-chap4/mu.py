#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=12
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
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    lZ = 1*np.log(np.sum(np.power(w,ly)))
    return -1/(ly*beta)*lZ


############################
# Déclaration matrices
J = 1 ; 

T = 1

muSpace = np.logspace(-3,0,20)
lSpace = np.arange(8,18)

energies = np.empty((len(lSpace),len(muSpace)))
casimir = np.empty((len(lSpace),len(muSpace)))

for i,mu in enumerate(muSpace) :
    for j,L1 in enumerate(lSpace):
        print(i,j)
        energies[j][i] = intDiag(L1,1/T,mu,J) 

for i,mu in enumerate(muSpace) :
    casimir[:,i] = -np.gradient(energies[:,i])


# Plot en fonction du potentiel, à L = 10
for i,L in enumerate(lSpace) :
    if L not in [10,15] : 
        continue
    plt.plot(muSpace,casimir[i,:],label="L="+str(L))

plt.xlabel('$\mu$')
plt.ylabel('$-\\frac{\partial \Omega(\\beta,L,h)}{\partial L}$')
plt.legend()
plt.tight_layout()
plt.savefig('casimir-mu-t'+str(T)+'.pdf')
plt.show()

