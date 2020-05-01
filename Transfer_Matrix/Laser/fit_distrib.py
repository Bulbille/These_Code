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
ncol = len(colors)
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
### Définition des matrices de magnetisation des différents modèles
def mat(L):
    d=np.zeros((L,L))
    for y,i in enumerate(d):
        d[y][y] = abs(y-L/2)
    return d
### Fonction qui calcule directement l'énergie libre
def distrib(ly,beta,h,j,mod):
    if mod == 'SOS' :
        tm = Trans(beta,h,j,ly)
    else :
        tm = TransGauss(beta,h,j,ly)
    w,v = LA.eigh(tm) 
    Z = np.sum(np.power(w,ly))
    ph = np.zeros(ly)
    for l in np.arange(ly):
        ph[l] = np.sum(np.power(w,ly)*np.power(v[l,:],2))
    return ph/Z

alpha   =sp.ai_zeros(3)[1]

def p(h,zeros,mean,l):
    nominateur = np.power(sp.airy(abs((h-mean)/l)+zeros[0])[0],2)
    denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(x/l+zeros[0])[0],2),0,np.inf)[0]
    return nominateur/denominateur

def pow(x,b,a):
    return b*np.power(x,a)

############################
# Déclaration matrices
J = 1 
lMax = 200
lSpace = np.arange(lMax)
fr = 30
xfit = lSpace[int(lMax/2)-fr:int(lMax/2)+fr+1]


tSpace = np.linspace(1,10,12)
bSpace = np.logspace(-1,0,1)
longueurs = np.empty((len(bSpace),len(tSpace)))

for i,t in enumerate(tSpace):
    for j,B in enumerate(bSpace):
        ph = distrib(lMax,1/t,B,J,'Gauss')
        phfit = ph[int(lMax/2)-fr:int(lMax/2)+fr+1]
        popt,pcov= curve_fit(lambda lSpace,l : np.log(p(lSpace,alpha,lMax/2,l)), xfit,np.log(phfit),p0=[3])
        longueurs[j][i] = abs(popt[0])
        print(t,B,popt,pcov)

#Plot en fonciton du champ magnétique à température constante
#Bexp = np.empty(len(tSpace))
#for i,t in enumerate(tSpace):
##    plt.plot(bSpace,longueurs[:,i],label="$T="+str(t)+'$')
#    popt,pcov = curve_fit(pow,bSpace,longueurs[:,i],p0=[1,-1/3])
##    plt.plot(bSpace,pow(bSpace,*popt),'+')
#    Bexp[i] = popt[1]
#    print('T',t,popt)
#plt.plot(tSpace,Bexp)
#plt.xlabel('$T$')
#plt.ylabel('Exposant de B')
#plt.legend()
#plt.savefig('exposantB.pdf')
#plt.show()
#exit()
#Plot en fonciton de la température à champ magnétique constant
Texp = np.empty(len(bSpace))
for i,b in enumerate(bSpace):
    plt.plot(tSpace,longueurs[i,:],label="$B="+str(b)+'$')
    popt,pcov = curve_fit(pow,tSpace,longueurs[i,:],p0=[1,2/3])
    plt.plot(tSpace,pow(tSpace,*popt),'+')
    print('B',b,popt)
    Texp[i] = popt[1]
plt.plot(bSpace,Texp)
plt.xlabel('$B$')
plt.xscale('log')
plt.ylabel('Exposant de T')
plt.legend()
plt.savefig('exposantT.pdf')
plt.show()

