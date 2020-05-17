#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=10
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from numpy import linalg as LA
import math as mt
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
### Définition des matrices de transfert des différents modèles
def TransSOS(kbt,champ,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(abs(y1)+abs(y2))/2 +inter*abs(y1-y2)) )
    return d
def TransPOP(kbt,champ,inter,L,mu):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*(-mu*(y1+y2)/2+ champ*(abs(y1)+abs(y2))/2
                +inter*abs(y1-y2))+(mt.log(mt.factorial(y1))+mt.log(mt.factorial(y2)))/2 )
    return d
### Définition des matrices de magnetisation des différents modèles
def Mat(L):
    d=np.zeros((L+1,L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(y)
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j,mod,mu):
    matphi = Mat(LY)
    matphi2 = np.square(matphi)
    if mod == 'SOS':
        tm = TransSOS(beta,h,j,ly)
    elif mod == 'POP' :
        tm = TransPOP(beta,h,j,ly,mu)
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
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1;

mumax = 5
muN = 10
####### Données #####
#### Calcul de la magnétisation via TM


hauteur = plt.subplot(311)
correlation = plt.subplot(312)
energie = plt.subplot(313)

mubar = np.linspace(0,mumax,muN)
hbar = np.empty(muN)
fbar = np.empty(muN)
cbar = np.empty(muN)
for s,string in enumerate(['POP']):
    for j,LY in enumerate(np.linspace(50,100,2).astype(int)) :
        print(LY,mubar)
        for k,mu in enumerate(mubar) :
            res = intDiag(LY,BETA,0,J,string,mu)
            fbar[k] = res[0]
            hbar[k] = res[1]
            cbar[k] = (res[2]-hbar[k]**2)**0.5

        hauteur.plot(mubar,hbar,color=colors[j],label=str(LY)+' '+str(string))
        hauteur.plot(mubar,np.exp(mubar),'+',color=colors[j])
        correlation.plot(mubar,cbar,color=colors[j])
        energie.plot(mubar,fbar,color=colors[j])

    
correlation.legend()
correlation.set_xlabel('$J$')
correlation.set_ylabel('Correlation length')

hauteur.legend()
hauteur.set_xlabel('$J$')
hauteur.set_ylabel('Height $\\bar{H}$')

energie.legend()
energie.set_xlabel('$J$')
energie.set_ylabel('$E=\sum |h_i-h_{i+1}|$')

plt.savefig('semiinf.pdf')
plt.show()
