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
    return a*np.power(x,-b)
### Définition des matrices de transfert des différents modèles
def Trans(kbt,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
           # d[y1][y2] = np.exp(-kbt*inter*abs(y1-y2)) 
            d[y1][y2] = np.exp(-kbt*inter*np.power(y1-y2,2)) 
    return d
### Fonction qui calcule directement l'énergie libre

def intDiag(ly,beta,j):
    tm = Trans(beta,j,ly)
    w,v = LA.eigh(tm) 
    for i in np.arange(3):
        plt.plot(np.arange(ly),v[:,i]**2)
        print(np.sum(v[:,i]**2))
#    plt.plot(np.arange(ly),pow(2/ly,0.5)*np.sin(np.pi*np.arange(ly)/ly),'+')
    return w[-1],-np.log(w[-1])/beta,-1/np.log(w[-2]/w[-1])


############################
# Déclaration matrices
J = 1 ; 

lMin =  101; lMax = 500
lSpace = np.arange(lMin,lMax)[::1000]
tSpace = np.linspace(2,10,1)

def fit(x,a,b):
    return a*np.power(x,-b)

for i,t in enumerate(tSpace) :
    ene = np.empty(len(lSpace))
    for j,L1 in enumerate(lSpace):
        res = intDiag(L1,1/t,J) 
        ene[j] = res[1]
        print(t,L1)
#    plt.plot(lSpace,ene,label=str(t))
#    plt.plot(lSpace,np.sinh(J/t)/(np.cosh(J/t)-np.cos(np.pi/lSpace)),'+')

plt.legend()

plt.show()
