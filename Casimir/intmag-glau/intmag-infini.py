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
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
def Trans2(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*(-champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def mat(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = y
    return d
def mat2(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = abs(y-L/2)
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    matphi = mat(ly)
    w,v = LA.eigh(tm) 
    try :
        phi = 0
        for i in np.arange(ly):
            phi+=np.dot(np.dot(v[:,i],matphi),v[:,i])*np.power(w[i],ly)
        phi /= np.sum(np.power(w,ly))
        return phi
    except :
        return np.dot(np.dot(v[:,-1],matphi),v[:,-1])
def intDiag2(ly,beta,h,j):
    tm = Trans2(beta,h,j,ly)
    matphi = mat2(ly)
    w,v = LA.eigh(tm) 
    return np.dot(np.dot(v[:,-1],matphi),v[:,-1])
    try:
        phi = 0
        for i in np.arange(ly) :
            phi+=np.dot(np.dot(v[:,i],matphi) , v[:,i])*np.power(w[i],ly)
        phi /= np.sum(np.power(w,ly))
        return phi
    except :
        return np.dot(np.dot(v[:,-1],matphi),v[:,-1])


############################
# Déclaration matrices
J = 1 ;  BETA  = 1
LY = 20
mumin = 1e-2 ; mumax = 1; muN = 100
muspace = np.linspace(mumin,mumax,muN)

mag = np.empty(muN)

for i,mu in enumerate(muspace):
    res= intDiag2(LY,BETA,mu,J)
    mag[i] = res

plt.plot(muspace,mag,label="Matrice de transfert")

data = np.loadtxt(sys.argv[1])

#integrale = np.zeros(len(data[:,1]))
#for i in np.arange(len(data[:,1])):
#    if i == 0 or i == len(data[:,1])-1 :
#        integrale[i] 


plt.plot(data[:,0],data[:,1],'+',label="Simulation numérique")

plt.xlabel('$\mu$')
plt.ylabel('$\langle \sum_i |h_i-\\frac{L}{2}| \\rangle $')
plt.legend()
plt.savefig('simu-tm-negstagged-l'+str(LY)+'.pdf')
plt.show()
