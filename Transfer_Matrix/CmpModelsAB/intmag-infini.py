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
def Trans2(kbt,champ,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*(-champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*abs(y1-y2)) )
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
    return -np.log(max(w))/beta
    lZ = 1*np.log(np.sum(np.power(w,ly)))
    return -1/(ly*beta)*lZ
def intDiag2(ly,beta,h,j):
    tm = Trans2(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    return -np.log(max(w))/beta
    lZ = 1*np.log(np.sum(np.power(w,ly)))
    return -1/(ly*beta)*lZ


############################
# Déclaration matrices
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1.;

LY = 20

hSpace = np.linspace(0,1,50)

def fitexp(x,a,b):
    return a*np.exp(-x/b)

data = np.empty(len(hSpace))
for i,h in enumerate(hSpace) :
    data[i] = intDiag(LY,BETA,h,J) 
plt.plot(hSpace,data,'+',label='$V(h_i) = h_i$')

for i,h in enumerate(hSpace) :
    data[i] = intDiag2(LY,BETA,h,J) 
plt.plot(hSpace,data+hSpace*LY/2,label='$V(h_i) = -|h_i-\\frac{L}{2}|$')

plt.xlabel('$\mu$')
plt.ylabel('$\Omega(\mu)-\Omega(\infty)$')

plt.legend()
plt.savefig('free-ene-potentiels.pdf')
plt.legend()
plt.show()
