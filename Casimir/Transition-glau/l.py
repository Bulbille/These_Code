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
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    lZ = 1*np.log(np.sum(np.power(w,ly)))
    return -1/(ly*beta)*lZ


############################
# Déclaration matrices
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1/TC;

data = np.loadtxt(sys.argv[1])
plt.plot(data[:,0],data[:,1],label="Simulation Glauber")

lSpace = data[:,0].astype(int)
print(lSpace)

h = 0.1
casimir = np.empty(len(lSpace))

for j,L1 in enumerate(lSpace):
    f1 = intDiag(L1,BETA,h,J) 
    casimir[j] = -f1
deriv = np.gradient(casimir)
plt.plot(lSpace,deriv,label="Matrice de Transfert")


plt.xlabel('$L$')
plt.ylabel('$f_c(\\beta,L,h)$')

plt.legend()
plt.savefig('casimir-temperature-zoom.pdf')
plt.show()
