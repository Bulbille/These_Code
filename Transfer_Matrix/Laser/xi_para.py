#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=12
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from numpy import linalg as LA
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate
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
            d[y1][y2] = np.exp(-kbt*( champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*abs(y1-y2)) )
    return d
def TransGauss(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*(y1-y2)**2) )
    return d
### Fonction qui calcule directement la longueur de corrélation
def Diag(ly,beta,h,j,mod):
    if mod == 'SOS' :
        tm = Trans(beta,h,j,ly)
    else :
        tm = TransGauss(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    return -1/np.log(w[-2]/w[-1])

def pow(x,b,a):
    return b*np.power(x,a)

############################
# Déclaration matrices
J = 1 
lMax = 200
lSpace = np.arange(lMax)

tSpace = np.linspace(4,10,10)
bSpace = np.logspace(-0,-1,2)
longueurs = np.empty((len(bSpace),len(tSpace)))

for i,t in enumerate(tSpace):
    for j,B in enumerate(bSpace):
        longueurs[j][i]= Diag(lMax,1/t,B,J,'Gauss')
        print(i,j,longueurs[j][i]**3/t)

##Plot en fonciton du champ magnétique à température constante
#Bexp = np.empty(len(tSpace))
#for i,t in enumerate(tSpace):
##    plt.plot(bSpace,longueurs[:,i],label="$T="+str(t)+'$')
#    popt,pcov = curve_fit(pow,bSpace,longueurs[:,i],p0=[1,-2/3])
##    plt.plot(bSpace,pow(bSpace,*popt),'+')
#    Bexp[i] = popt[1]
#    print('T',t,popt)
#plt.plot(tSpace,Bexp)
#plt.xlabel('$T$')
#plt.ylabel('Exposant de B')
#plt.legend()
#plt.savefig('exposantB-para.pdf')
#plt.show()
#Plot en fonciton de la température à champ magnétique constant
Texp = np.empty(len(bSpace))
for i,b in enumerate(bSpace):
    plt.plot(tSpace,longueurs[i,:],label="$B="+str(b)+'$')
    popt,pcov = curve_fit(pow,tSpace,longueurs[i,:],p0=[1,1/3])
    plt.plot(tSpace,pow(tSpace,*popt),'+')
    print('B',b,popt)
    Texp[i] = popt[1]
#plt.plot(bSpace,Texp)
plt.xlabel('$B$')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Exposant de T')
plt.legend()
plt.savefig('exposantT-para.pdf')
plt.show()

