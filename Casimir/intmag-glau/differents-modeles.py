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
def TransAbs(kbt,champ,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(abs(L/2-y1)+abs(L/2-y2))/2 +inter*abs(y1-y2)) )
    return d

def TransSOS(kbt,champ,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def MatSOS(L):
    d=np.zeros((L+1,L+1))
    for y,i in enumerate(d):
        d[y][y] = y
    return d
def MatAbs(L):
    d=np.zeros((L+1,L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(L/2-y)
    return d
### Fonction qui calcule directement l'énergie libre
def DiagSOS(ly,beta,h,j):
    matphi = MatSOS(LY)
    tm = TransSOS(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    return [np.dot( np.dot(v[:,np.argmax(w)],matphi) , v[:,np.argmax(w)]) , -1/beta*np.log(np.sum(np.power(w,ly)))]

def DiagAbs(ly,beta,h,j):
    matphi = MatAbs(LY)
    tm = TransAbs(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    return [np.dot( np.dot(v[:,np.argmax(w)],matphi) , v[:,np.argmax(w)]) , -1/beta*np.log(np.sum(np.power(w,ly)))]
############################
# Déclaration matrices
J = 1 ; BETA=1;

Bspace = np.linspace(1e-2,1,50)
####### Données #####
#### Calcul de la magnétisation via TM

hsos = np.empty(np.size(Bspace))
habs = np.empty(np.size(Bspace))

hauteur = plt.subplot(111)

LY = 40
for i,Bmax in enumerate(Bspace):
    res = DiagSOS(LY,BETA,Bmax,J);
    hsos[i] = res[1];
    res = DiagAbs(LY,BETA,Bmax,J);
    habs[i] = res[1]; 

hauteur.plot(Bspace,hsos,label="SOS")
hauteur.plot(Bspace,habs,label="Abs")

plt.show()
