#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=10
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import brentq
from numpy import linalg as LA
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
### Définition des matrices de transfert des différents modèles
def Trans(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*inter*abs(y1-y2))
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) ; 
    lZ = 1*np.log(max(w))
    return -1/beta*lZ

#def equaDiag(ly,beta,j):
#    betexp = lambda beta : np.exp(-beta)
#    theta = lambda eig,h,r : 1/np.tan(h*eig) - r*np.sin(eig)/(1-r*np.cos(eig))
#    ret:

############################
# Déclaration matrices
TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1;

axfree = plt.subplot(111)

def dean(j,ly):
    return - np.log(np.sinh(j)/(np.cosh(j)-np.cos(np.pi/(2*ly+2))))

Jn = 30
Jspace = np.round(np.linspace(0.3,3,Jn),2)
for nl,LY in enumerate(int(np.linspace(10,70,4))):
    FreTM = np.empty([np.size(Jspace),2])
    for nb,J in enumerate(Jspace):
        FreTM[nb] = [J,intDiag(LY,BETA,0,J)]

    if nl == 0 :
        axfree.plot(FreTM[:,0],FreTM[:,1],color=colors[nl%nb_col],label="Numerical : $LY="+str(LY)+'$')
        axfree.plot(FreTM[:,0],dean(FreTM[:,0],LY),'+',color=colors[nl%nb_col],label="Analytic approximation")
    else :
        axfree.plot(FreTM[:,0],FreTM[:,1],color=colors[nl%nb_col],label="$LY="+str(LY)+'$')
        axfree.plot(FreTM[:,0],dean(FreTM[:,0],LY),'+',color=colors[nl%nb_col])

axfree.legend(loc='upper right')
axfree.set_xlabel('$J$')
axfree.set_ylabel('$F(L_Y,J)$')

plt.savefig('null_deanJ.pdf')
plt.show()
