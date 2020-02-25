#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fontsize=12
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import scipy as sy
from scipy.optimize import curve_fit
import numpy as np
from numpy import linalg as LA
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab','palegreen','lightseagreen','darkcyan','royalblue','navy','plum','coral','orange','olive','g','c','dodgerblue','violet','fuchsia','crimson','peru','goldenrod','forestgreen','aquamarine','aqua','blueviolet','pink','tomato','linen','antiquewhite','yellow','greenyellow','lime','cyan','indigo','darkmagenta']


################################## 
### Déclaration fonctions ########
##################################

### Définition des matrices de transfert des différents modèles
def TransA(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*champ*((L-y1)+(L-y2))/2 -kbt*inter*abs(y1-y2) )
    return d
def TransB(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(abs(L-y1)+abs(L-y2))/2 +inter*abs(y1-y2)) )
    return d
def TransC(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( -champ*(abs(L-y1)+abs(L-y2))/2 +inter*abs(y1-y2)) )
    return d
### Définition des matrices de magnetisation des différents modèles
def matA(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = (L-y)
    return d
def matB(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = abs(L-y)
    return d
def matC(L):
    d=np.zeros((2*L+1,2*L+1))
    for y,i in enumerate(d):
        d[y][y] = -abs(L-y)
    return d
### Fonction qui intègre de 0 à Bmax.
# numero : 0 pour modèle A, 1 pour modèle B, 2 pour modèle C
def intMag(lx,ly,hmin,hmax,hn,beta,j,numero):
    omega = 0
    if numero == 0:
        matphi = matA(ly)
    elif numero == 1:
        matphi = matB(ly)
    else :
        matphi = matC(ly)

    for k,h in enumerate(np.linspace(hmin,hmax,hn)):
        if numero == 0:
            tm = TransA(beta,h,j,ly)
        elif numero == 1:
            tm = TransB(beta,h,j,ly)
        else :
            tm = TransC(beta,h,j,ly)

        w,v = LA.eigh(tm)
        # Limite LX → infty
        # m = <L0|D|L0>
        phi = np.dot( np.dot(v[:,np.argmax(w)],matphi) , v[:,np.argmax(w)])

        if(k == 0 or k == hn-1):
            omega += phi
        elif (k%2 == 0 ):
            omega += 2*phi
        else:
            omega += 4*phi
    omega /= (3*hn)/(hmax-hmin)
    return omega
### Fonction qui calcule directement l'énergie libre
def intDiag(lx,ly,beta,h,j,numero):
    if numero == 0:
        tm = TransA(beta,h,j,ly)
    elif numero == 1:
        tm = TransB(beta,h,j,ly)
    else :
        tm = TransC(beta,h,j,ly)
    w,v = LA.eigh(tm) ; lZ = 1*np.log(max(w))
    return -1/beta*lZ
############################
# Déclaration matrices
BETA = 1; J = 1; 
LX = 70
LY = 30

Hspace = 10
mags    = np.linspace(0,0.2,Hspace)

IntA = np.zeros(Hspace)
IntB = np.zeros(Hspace)
IntC = np.zeros(Hspace)

Hmin = 0
Hn   = 100

FreA = np.empty(Hspace)
FreB = np.empty(Hspace)
FreC = np.empty(Hspace)
FreA[:] = intDiag(LX,LY,BETA,Hmin,J,0)
FreB[:] = intDiag(LX,LY,BETA,Hmin,J,1)
FreC[:] = intDiag(LX,LY,BETA,Hmin,J,2)


for nb,Hmax in enumerate(mags):
    print nb

#    IntA[nb:] -= intMag(LX,LY,Hmin,Hmax,Hn,BETA,J,0);
#    FreA[nb]  -= intDiag(LX,LY,BETA,Hmax,J,0)

    IntB[nb:] -= intMag(LX,LY,Hmin,Hmax,Hn,BETA,J,1);
    FreB[nb]  -= intDiag(LX,LY,BETA,Hmax,J,1)

    IntC[nb:] -= intMag(LX,LY,Hmin,Hmax,Hn,BETA,J,2);
    FreC[nb]  -= intDiag(LX,LY,BETA,Hmax,J,2)

    Hmin = Hmax

#plt.plot(mags,IntA,'+',color='blue',label='Model A : magnetisation')
#plt.plot(mags,FreA,'-.',color='blue',label='Model A : diagonalisation')

plt.plot(mags,IntB,'+',color='red',label='Model B : magnetisation')
plt.plot(mags,FreB,':',color='red',label='Model B : diagonalisation')

plt.plot(mags,IntC,'+',color='green',label='Model C : magnetisation')
plt.plot(mags,FreC,':',color='green',label='Model C : diagonalisation')

def fit_lin(x,a):
    return -a*x

#popt,pcov = curve_fit(fit_lin,mags,IntA)
#print popt
#popt,pcov = curve_fit(fit_lin,mags,FreA)
#print popt

plt.legend(loc='upper left')
plt.ylabel('Free energy')
plt.xlabel('$B^\\ast$')
plt.savefig('comparison.pdf')
plt.show()
plt.close()


