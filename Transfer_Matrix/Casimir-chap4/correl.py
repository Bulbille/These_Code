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

def fitexp(x,a,l):
    return a*np.exp(-x/l)

### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    Z = np.sum(np.power(w,ly))

    matphi = mat(ly)
    correlthermo = -1/np.log(w[-2]/w[-1])
    return [correlthermo]


############################
# Déclaration matrices
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1/TC;

LY = 200

muSpace = np.logspace(-4,1,120)
tSpace = np.linspace(0.4,2,3)
casimir = np.empty(len(muSpace))

correl = np.empty((len(muSpace),1))

def fitexp(x,a,l):
    return a*np.power(x,l)

for i,t in enumerate(tSpace) :
    for j,mu in enumerate(muSpace):
        f1 = intDiag(LY,1/(t*TC),mu,J) 
        correl[j] = f1
    plt.plot(muSpace,correl[:,0],label='T='+str(t)[:3])
    popt,pcov = curve_fit(fitexp,muSpace[40:80],correl[40:80,0],p0=[1,-1])
    print(popt)
#    plt.plot(muSpace,fitexp(muSpace,*popt))

plt.xlabel('$\mu$')
plt.ylabel('$\\xi=- \\frac{1}{- \ln(\\frac{\lambda_1}{\lambda_0})}$')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.savefig('longueur-correl.pdf')
plt.show()

