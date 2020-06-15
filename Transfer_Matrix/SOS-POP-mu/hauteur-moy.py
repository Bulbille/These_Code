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
import math as mt
import scipy.special as sp
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
linestyles = ['-','-.',':','--']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
### Définition des matrices de transfert des différents modèles
def TransSOS(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*(y1-y2)**2)/2 )
    return d


### Définition des matrices de magnetisation des différents modèles
def Mat(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = y
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    matphi = Mat(ly)
    tm = TransSOS(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    phi =  np.dot(np.dot(v[:,-1],matphi),v[:,-1])
    return phi,-np.log(w[-1])/beta

alpha =sp.ai_zeros(4)[1][0]
def hmoy(sig,beta,mu) :
    return 2/3*alpha/np.power(2*sig*beta**2*mu,1/3)

############################
# Déclaration matrices
J = 2 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1;

Bn = 12
Bspace = np.linspace(0,3,Bn)
Bspace = [0.0001]
####### Données #####
#### Calcul de la magnétisation via TM

hbar = np.empty(np.size(Bspace))

hauteur = plt.subplot(111)

LY = 100

for i,Bmax in enumerate(Bspace) :
    res = intDiag(LY,BETA,Bmax,J)
    hbar[i] = res[0]
    print(Bmax,res[0],-hmoy(J,BETA,Bmax))

hauteur.plot(Bspace,hbar,label="Transfer matrix")
#hauteur.plot(Bspace,-hmoy(J,BETA,Bspace),label="Continuous model")

hauteur.legend()
hauteur.set_xlabel('$\mu$')
#hauteur.set_xscale('log')
##hauteur.set_yscale('log')
hauteur.set_ylabel('$\langle h \\rangle $')

plt.savefig('hauteur-tm-sos.pdf')
plt.show()
