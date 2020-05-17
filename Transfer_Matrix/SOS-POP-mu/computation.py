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
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
linestyles = ['-','-.',':','--']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
### Définition des matrices de transfert des différents modèles
def TransSOS(kbt,champ,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
def lnfact(n):
    return mt.log(mt.factorial(n))

def TransPOP(kbt,mu,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2]= np.exp(-kbt*(-mu*(y1+y2)/2+inter*abs(y1-y2))-(lnfact(y1)+lnfact(y2))/2)
    return d

### Définition des matrices de magnetisation des différents modèles
def Mat(L):
    d=np.zeros((L+1,L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(y)
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j,mod):
    matphi = Mat(LY)
    matphi2 = np.square(matphi)
    if mod == 'SOS':
        tm = TransSOS(beta,h,j,ly)
    elif mod == 'POP' :
        tm = TransPOP(beta,h,j,ly)   
    w,v = LA.eigh(tm) 
    phi =  np.dot(np.dot(v[:,-1],matphi),v[:,-1])
    phi2 =  np.dot(np.dot(v[:,-1],matphi2),v[:,-1])
    sigma = (phi2-phi**2)**0.5

    return [phi,sigma,-1/beta*np.log(max(w))]

############################
# Déclaration matrices
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1;

Bn = 60
Bspace = np.logspace(-1,0.77,Bn)
####### Données #####
#### Calcul de la magnétisation via TM

hbar = np.empty(np.size(Bspace))
fbar = np.empty(np.size(Bspace))
cbar = np.empty(np.size(Bspace))

hauteur = plt.subplot(111)

def exp(x,a,l):
    return a*np.exp(-x/l)

for s,string in enumerate(['POP']):
    for l,LY in enumerate([25,75,200][::-1]) :
        print(s,LY)
        for i,Bmax in enumerate(Bspace) :
            res = intDiag(LY,BETA,Bmax,J,string)
            hbar[i] = res[0]
            cbar[i] = res[1]
            fbar[i] = res[2]

        hauteur.plot(Bspace,hbar,color=colors[l],label="L="+str(LY))

hauteur.plot(Bspace[:80],np.exp(BETA*Bspace[:80]),'+',color="black",label="$e^{\\beta \mu}$")

hauteur.legend()
hauteur.set_xlabel('$\mu$')
hauteur.set_xscale('log')
hauteur.set_yscale('log')
hauteur.set_ylabel('$\langle h \\rangle $')

plt.savefig('hauteur-tm-pop.pdf')
plt.show()
