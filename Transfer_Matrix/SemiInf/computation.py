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
### Définition des matrices de transfert des différents modèles
def V(y) :
    if y==0 :
        return 1e8
    return y

def Trans(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(V(y1)+V(y2))/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def Mat(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = abs(y)
    return d
### Fonction qui calcule directement l'énergie libre
def p(ly,beta,h,j):
    matphi = Mat(LY)
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    return np.dot(np.dot(v[:,-1],matphi),v[:,-1])
    return np.power(v[:,-1],2)

import scipy.special as sp
alpha =sp.ai_zeros(1)[0][0]
print(alpha)
def hmoy(sig,beta,mu) :
    return 2/3*alpha/np.power(2*sig*beta**2*mu,1/3)



############################
# Déclaration matrices
J = 1 ; BETA = 1;
LY = 100
bs = np.linspace(0.01,1,10)
ana = np.empty(len(bs))
gsos = np.empty(len(bs))

for i,B in enumerate(bs) :
    gsos[i] = p(LY,BETA,B,J)
    ana[i] = hmoy(J/2,BETA,B)
plt.plot(bs,-ana,label="Ana")
plt.plot(bs,gsos,label="GSOS")
plt.legend()
plt.show()
