#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=14
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
import sys
from numpy import linalg as LA
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab','palegreen','lightseagreen','darkcyan','royalblue','navy','plum','coral','orange','olive','g','c','dodgerblue','violet','fuchsia','crimson','peru','goldenrod','forestgreen','aquamarine','aqua','blueviolet','pink','tomato','linen','antiquewhite','yellow','greenyellow','lime','cyan','indigo','darkmagenta']


################################## 
### Déclaration fonctions ########
##################################

### Définition des matrices de transfert des différents modèles
def TransA(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def matA(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = y
    return d
### Fonction qui intègre de 0 à Bmax.
# numero : 0 pour modèle A, 1 pour modèle B, 2 pour modèle C
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = TransA(beta,h,j,ly)
    matphi = matA(ly)
    w,v = LA.eigh(tm) 
    return np.power(v[:,-1],2)
############################
# Déclaration matrices
BETA = 1; J = 1; 

for LY in [10,80] :
    plt.plot(np.arange(LY),intDiag(LY,BETA,0,J),label="Transfer matrix, $L="+str(LY)+'$')

data=np.loadtxt(sys.argv[1])
plt.plot(data[:,0],data[:,1]/np.sum(data[:,1]),label="Kawasaki, $L=10$")
#data=np.loadtxt(sys.argv[2])
#plt.plot(data[:,0],data[:,1]/np.sum(data[:,1]),label="Kawasaki, $L=10$")

plt.xlabel('$h$')
plt.ylabel('$p(h)$')
plt.legend(loc='upper right')
plt.savefig('pdf-kaw.pdf')
plt.show()
