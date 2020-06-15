#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=12
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
    return a*np.power(x,-b)
### Définition des matrices de transfert des différents modèles
def Trans(kbt,champ,inter,L):
    d=np.zeros((L+1,L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
#            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*pow(y1-y2,2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def mat(L):
    d=np.zeros((L+1,L+1))
    for y,i in enumerate(d):
        d[y][y] = y
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    plt.plot(np.arange(ly),v[:,-1])
    return -np.log(w[-1])/beta


############################
# Déclaration matrices
J = 1 ; 

lMin = 10; lMax = 30
lSpace = np.arange(lMin,lMax)
tSpace = np.linspace(3.6,5.3,1)

h = 0.0

energies = np.empty((len(lSpace),len(tSpace)))
casimir = np.empty((len(lSpace),len(tSpace)))

for i,t in enumerate(tSpace) :
    for j,L1 in enumerate(lSpace):
        energies[j][i] = intDiag(L1,1/t,h,J) 
        plt.scatter(t**2,energies[j][i])

#for i,t in enumerate(tSpace) :
#    casimir[:,i] = -np.gradient(energies[:,i])
#    casimir[:,i] = energies[:,i]
#
# Plot en fonction de la distance, à température constante
# Ne pas montrer premier et dernier élément des dérivées, moins bonne qualité
#for i,t in enumerate(tSpace) :
##     if i % 4 != 0 :
##         continue
#    plt.plot(lSpace[1:-1],casimir[1:-1,i],label="T="+str(t)[:4])
#    popt,pcoc = curve_fit(linear,lSpace[1:-1],casimir[1:-1,i],p0=[-1,3])
#    print(popt)
#    plt.plot(lSpace[1:-1],linear(lSpace[1:-1],*popt),'+')


plt.xlabel('$L$')
plt.ylabel('$-\\frac{\partial \Omega(\\beta,L,h)}{\partial L}$')
plt.legend()
plt.tight_layout()
plt.savefig('casimir-distance-mu'+str(h)+'.pdf')
plt.show()
exit()

# Plot en fonction de la température, à distance constante
for i,L in enumerate(lSpace) :
    if L not in [5,10,15,20] :
        continue
    plt.plot(tSpace,casimir[i,:],label="L="+str(L))

plt.xlabel('$T$')
plt.ylabel('$-\\frac{\partial \Omega(\\beta,L,h)}{\partial L}$')
plt.legend()
plt.tight_layout()
plt.savefig('casimir-temperature-mu'+str(h)+'.pdf')
plt.show()

