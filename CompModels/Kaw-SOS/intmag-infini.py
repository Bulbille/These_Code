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
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def mat(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = y
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    Z = np.sum(np.power(w,200))
    matphi = mat(LY)
    phi = 0
    for i,l in enumerate(w):
        phi+= np.dot(np.dot(v[:,i],matphi),v[:,i])*np.power(l,200)

#    return phi/Z
    return np.dot(np.dot(v[:,-1],matphi),v[:,-1])


############################
# Déclaration matrices
J = 1 ;BETA = 1;

LY = 200

hSpace = np.linspace(0.01,1,8)


data = np.empty(len(hSpace))
for i,h in enumerate(hSpace) :
    data[i] = intDiag(LY,BETA,h,J) 
plt.plot(hSpace,data,'+-',label='Matrice de transfert')
print(data)

data = np.loadtxt(sys.argv[1])
plt.plot(data[:,0],data[:,1],label="Glauber")
print(data[:,1])
plt.xlabel('$\mu$')

plt.legend()
plt.show()
