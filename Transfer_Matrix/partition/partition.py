#!/usr/bin/python
# -*- coding: utf8
import numpy as np

import matplotlib.pyplot as plt
from numpy import linalg as LA
fontsize=10
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

################################## 
### Déclaration fonctions ########
##################################

#Coefficient multinomiaux 
#->sert pour calculer le poids statistique d'une configuration 
from scipy.special import binom
def multinomial(params):
    if len(params) == 1:
        return 1
    return binom(sum(params), params[-1]) * multinomial(params[:-1])

##################################
def Trans(kbt,champ,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm)
    return np.sum(np.power(w,ly))

##################################

TC  = 2/np.log(1.+np.power(2,0.5)); BETA = 1/TC ; J = 1; 

LY = 1
Zc = 0
Zg = 0
for M in np.arange(0,2*LY+1) :
    for h1 in np.arange(0,M+1):
        for h2 in np.arange(0,M+1-h1) :
            for h3 in np.arange(0,M+1-h1-h2):
                if h1+h2+h3 == LY:
                    Zc += np.exp(-BETA*(abs(h1-h2)+abs(h2-h3)+abs(h3-h1)))

print(Zc)

yrange = np.arange(0,2*LY+1)
for h1 in yrange:
    for h2 in yrange :
        for h3 in yrange:
            Zg += np.exp(-BETA*(abs(h1-h2)+abs(h2-h3)+abs(h3-h1)))
print(Zg)
