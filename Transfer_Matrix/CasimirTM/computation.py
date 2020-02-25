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
def Trans(kbt,champ,inter,L,numero):
    if numero == 'A':
        return TransA(kbt,champ,inter,L)
    elif numero =='B' :
        return TransB(kbt,champ,inter,L)
    else :
        return TransC(kbt,champ,inter,L)
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
def Mat(L,numero):
    if numero == 'A':
        return matA(L)
    elif numero =='B' :
        return matB(L)
    else:
        return matC(L)
### Fonction qui intègre de 0 à Bmax.
# numero : 0 pour modèle A, 1 pour modèle B, 2 pour modèle C
def intMag(lx,ly,hmin,hmax,hn,beta,j,numero):
    omega = 0
    matphi = Mat(ly,numero)
    for k,h in enumerate(np.linspace(hmin,hmax,hn)):
        tm = Trans(beta,h,j,ly)
        w,v = LA.eigh(tm)
        #Pas de limite
        Z = sum(w**lx)
        tmpow = LA.matrix_power(tm,LX)
        phi = np.trace(np.dot(matphi,tmpow))/Z

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
    tm = Trans(beta,h,j,ly,numero)
    w,v = LA.eigh(tm) ; 
    lZ = 1*np.log(max(w))
    return -1/beta*lZ

def asymptot(ly,beta,h,j,numero):
    if numero == 'A' :
        return -ly*h 
    elif numero == 'B' :
        return ly*0
    elif numero == 'C' :
        return -ly*h - 1/beta * np.log(1+np.exp(-2*beta*j*ly))
    else :
        return 0

############################
# Déclaration matrices
J = 1 ; TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 3/TC;
LX = 70
model = 'C'

Hn = 10
Hmax = 0.05
Hspace = np.round(np.linspace(0,Hmax,Hn),3)
Lspace = np.concatenate((np.arange(5,40),[70]))

####### Données #####
#### Calcul de la magnétisation via TM
def fitpow(x,l,a,b):
    return a*np.power(x,-l)+b
def fitpowexp(x,lam,alpha,a,b):
    return a*np.exp(-x/lam)*np.power(x,alpha)+b

lengthexp = np.empty(Hn)
lengthpow = np.empty(Hn)

for nb,Hmax in enumerate(Hspace):
    if nb%2 != 0 and nb != len(Hspace)-1:
        continue
    FreTM = np.empty([np.size(Lspace),2])
    for numero,LY in enumerate(Lspace):
        FreTM[numero] = [LY,intDiag(LX,LY,BETA,Hmax,J,model)]

    if model == 'B':
        FreTM[:,1] = FreTM[:,1]-FreTM[-1,1]
        Nmax = np.where(FreTM[:,1]>1e-14)[0][-1]
        plt.plot(FreTM[:Nmax,0],FreTM[:Nmax,1],label='$H = '+str(Hmax)+'$',color=colors[nb%nb_col])
    
        try :
            popw,pcov = curve_fit(fitpow,FreTM[:Nmax,0],FreTM[:Nmax,1],p0=[2,1,1])
            print popw
            plt.plot(FreTM[:Nmax,0],fitpow(FreTM[:Nmax,0],*popw),'o',color=colors[nb%nb_col])
        except:
            pass
        try :
            popw,pcov = curve_fit(fitpowexp,FreTM[:Nmax,0],FreTM[:Nmax,1],p0=[5,0,1,1])
            print popw
            plt.plot(FreTM[:Nmax,0],fitpowexp(FreTM[:Nmax,0],*popw),'x',color=colors[nb%nb_col])
        except :
            pass

    elif model == 'A' or model == 'C':
        FreTM[:,1] = FreTM[:,1]/FreTM[:,0]
        FreTM[:,1] = FreTM[:,1]-FreTM[-1,1]
        FreTM[:,1] *= -1
        popw,pcov = curve_fit(fitpow,FreTM[:,0],FreTM[:,1],p0=[2,0.1,0])
        print popw
        plt.plot(FreTM[:-1,0],FreTM[:-1,1],label='$H = '+str(Hmax)+'$',color=colors[nb%nb_col])
        plt.plot(FreTM[:-1,0],fitpow(FreTM[:-1,0],*popw),'x',color=colors[nb%nb_col])

plt.yscale('log')
#plt.xscale('log')
plt.legend(loc='lower right')
plt.xlabel('$L_Y$')
plt.ylabel('$F(L_Y,B)-F(\infty,B)$')

plt.title('Model '+model)
plt.savefig('free_energy'+model+'.pdf')
plt.show()
