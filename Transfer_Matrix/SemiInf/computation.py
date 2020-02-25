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
def Trans(kbt,champ,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def Mat(L):
    d=np.zeros((L+1,L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(y)
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    matphi = Mat(LY)
    matphi2 = np.square(matphi)
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 

    deb = np.dot(np.dot(v[:,-1],matphi),v[:,-2])*np.dot(np.dot(v[:,-1],matphi),v[:,-2])
    r = np.arange(512)
    r = deb*np.power(w[-2]/w[-1],r)
    chi = 0
    for i,c in enumerate(r):
        if i == 0 or i == len(r)-1:
            chi += c/r[0]
        elif i%2 == 0 :
            chi += 2*c/r[0]
        else :
            chi += 4*c/r[0]
    chi /=3
    return [-np.log(max(w))/beta , np.dot(np.dot(v[:,-1],matphi),v[:,-1]),np.dot(np.dot(v[:,-1],matphi2),v[:,-1]), chi] 

############################
# Déclaration matrices
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1/TC;

Bn = 50
Bspace = np.logspace(-6,1,Bn)
LY = 100
####### Données #####
#### Calcul de la magnétisation via TM

hbar = np.empty(np.size(Bspace))
fbar = np.empty(np.size(Bspace))
cbar = np.empty(np.size(Bspace))

hauteur = plt.subplot(111)
#correlation = plt.subplot(312)
#energie = plt.subplot(313)

for l,LY in enumerate([50,100,150,200]) :
    for i,Bmax in enumerate(Bspace) :
        print(l,i)
        res = intDiag(LY,BETA,Bmax,J)
        fbar[i] = res[0]
        hbar[i] = res[1]
        cbar[i] = (res[2]-hbar[i]**2)**0.5

    hauteur.plot(Bspace,hbar,color=colors[l],label="N="+str(LY))

    
#correlation.legend()
#correlation.set_xlabel('$\mu$')
#correlation.set_ylabel('Correlation length')

hauteur.legend()
hauteur.set_xlabel('$\mu$')
hauteur.set_xscale('log')
hauteur.set_ylabel('$\langle H \\rangle $')

#energie.legend()
#energie.set_xlabel('$\mu$')
#energie.set_ylabel('$E=\sum |h_i-h_{i+1}|$')

plt.savefig('semiinf.pdf')
plt.show()
