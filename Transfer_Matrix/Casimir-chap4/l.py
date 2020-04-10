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

lMin = 10; lMax = 80
lSpace = np.arange(lMin,lMax)

tSpace = np.linspace(0.2,1,10)
h = 0.1
casimir = np.empty(len(lSpace))

def fitpow(x,n,a):
    return a*np.power(x,n)
def fitlin(x,n,a):
    return a*x+n

for i,t in enumerate(tSpace) :
    for j,L1 in enumerate(lSpace):
        f1 = intDiag(L1,1/(t*TC),h,J) 
        casimir[j] = -f1
    deriv = np.gradient(casimir)
    plt.plot(lSpace,deriv,label=str(t),color=colors[i%len(colors)])

#    popt,pcov=curve_fit(fitpow,lSpace,deriv,p0=[-1,0])
#    print(t,popt)
#    plt.plot(lSpace,fitpow(lSpace,*popt),'+',color=colors[i%len(colors)])


plt.legend()
plt.show()
