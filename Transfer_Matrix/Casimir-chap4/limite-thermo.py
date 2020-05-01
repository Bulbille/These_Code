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
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(y1+y2)/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def mat(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] =y
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
#    matphi = mat(ly)
#    phi = 0
#    for i in np.arange(ly):
#        phi+=np.dot(np.dot(v[:,i],matphi) , v[:,i])*np.power(w[i],ly)
#    phi /= np.sum(np.power(w,ly))
#    return [max(w),phi,np.dot(np.dot(v[:,np.argmax(w)],matphi) , v[:,np.argmax(w)])]

    return max(w),-1/beta*np.log(max(w)), -1/(ly*beta)*np.log(np.sum(np.power(w,ly)))

J = 1 ; BETA = 1 ; 

LSpace = np.arange(5,50)
partition = np.empty((len(LSpace),4))
for i,L in enumerate(LSpace):
    res = intDiag(L,BETA,0.1,J)
    partition[i][3] = L
    partition[i][0] = res[0]
    partition[i][1] = res[1]
    partition[i][2] = res[2]

plt.figure(1)
freeEne = plt.subplot(111)

freeEne.plot(partition[:,3],partition[:,0],label='$\lambda_{0}$')
freeEne.set_xlabel('$L$')

plt.legend()
plt.savefig('freeene-lambda0-mu.pdf')
plt.show()

plt.figure(2)
thermo  = plt.subplot(111)

thermo.plot(partition[:,3],partition[:,1],label="$-\\frac{1}{\\beta} \ln(\lambda_0)$")
thermo.plot(partition[:,3],partition[:,2],label="$-\\frac{1}{L \\beta} \ln(\sum_\lambda \lambda^L)$")

np.savetxt("data",partition)

thermo.set_xlabel('$L$')
thermo.set_ylabel('$F(L,\\beta)$')

plt.legend()
plt.savefig('freeene-thermo-mu.pdf')
plt.show()
