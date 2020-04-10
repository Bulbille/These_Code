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
import time
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
def linear(x,a,b):
    return a*x+b
### Définition des matrices de transfert des différents modèles
def Trans(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def mat(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(y)
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    lZ = 1*np.log(max(w))
    return -1/beta*lZ


############################
# Déclaration matrices
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1/TC;

L1 = 10
L2 = 50

tN = 60; tMAX = 6; tMIN = 0.001
tSpace = np.linspace(tMIN,tMAX,tN) 
hN = 10; hMAX = 0.2; hMIN = 0
hSpace = np.linspace(hMIN,hMAX,hN) 
h = 0
casimir = np.empty(tN)
nu = 1
d= 1

def fitpow(x,n,a):
    return a*np.power(x,n)
def fitlin(x,n,a):
    return a*x+n

#for j,L1 in enumerate([5,10,15,20,25,50]):
for j,L1 in enumerate([16]):
    for i,t in enumerate(tSpace) :
#        f1 = intDiag(L1+1,1/(t*TC),h,J) 
        f1 = intDiag(L1,1/(t*TC),h,J) 
        casimir[i] = -f1
    deriv = np.gradient(casimir)

    plt.plot(tSpace,deriv,label=str(L1),color=colors[j])

#    popt,pcov=curve_fit(fitlin,tSpace[30:],casimir[30:],p0=[1,1])
#    plt.plot(tSpace,fitlin(tSpace,*popt),'+',color=colors[j])


plt.legend()
plt.show()
