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
        d[y][y] = abs(y)
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(ly,beta,h,j,lx):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    return w[-1],w[-2],-1/beta*np.log(max(w)), -1/(lx*beta)*np.log(np.sum(np.power(w,lx)))

J = 1 ; BETA = 1 ; 
LY = 100

LSpace = np.arange(5,150)
partition = np.empty((len(LSpace),5))
for i,L in enumerate(LSpace):
    print(L)
    res =  intDiag(LY,BETA,0,J,L)
    partition[i][3] = L
    partition[i][0] = res[0]
    partition[i][4] = res[1]
    partition[i][1] = res[2]
    partition[i][2] = res[3]

plt.figure(1)
freeEne = plt.subplot(111)

freeEne.plot(partition[:,3],partition[:,0],label='$\lambda_{0}$')
freeEne.plot(partition[:,3],partition[:,0],label='$\lambda_{1}$')
freeEne.plot(partition[:,3],0*partition[:,3]+ np.sinh(BETA*J)/(np.cosh(BETA*J)-1),'-.',color="black",label="$\\frac{\sinh(\\beta J)}{\cosh(\\beta J) -1}$")
freeEne.set_xlabel('$L$')

plt.legend()
plt.savefig('freeene-lambda0-libre.pdf')
plt.show()

plt.figure(2)
thermo  = plt.subplot(111)

thermo.plot(partition[:,3],partition[:,1],label="$F(\infty)= -\\frac{1}{\\beta} \ln(\lambda_0)$")
thermo.plot(partition[:,3],partition[:,2],label="$F(L') =-\\frac{1}{L\' \\beta} \ln(\sum_\lambda \lambda^{L\'})$")

np.savetxt("data",partition)

thermo.set_xlabel('$L\'$')
thermo.set_ylabel('$F(L\',\\beta)$')

plt.legend()
plt.savefig('freeene-thermo-libre.pdf')
plt.show()
