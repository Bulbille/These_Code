#!/usr/bin/python
# -*- coding: utf8

#########################################
# Modèle d'interface en histogramme
#########################################
import matplotlib.pyplot as plt
fontsize=9
labelsize=16
plt.rc("font",size=fontsize)
plt.rc("font",family='serif')

import numpy as np
import os #pour les chemins de fichiers
import sys #pour l'exit et l'argument
import re

from scipy.optimize import curve_fit
from numpy import linalg as LA
import scipy as sy
import scipy.special as sp
import scipy.integrate as integrate

colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
ncol = len(colors)


alpha   =sp.ai_zeros(3)[1][0]
alphap   =sp.ai_zeros(3)[0][0]

def p(h,zeros,mean,l):
    nominateur = np.power(sp.airy(abs((h-mean)/l)+zeros)[0],2)
    denominateur = 2*integrate.quad(lambda x : np.power(sp.airy(x/l+zeros)[0],2),0,np.inf)[0]
    return nominateur/denominateur
def TransGauss(kbt,champ,inter,L):
    d=np.zeros((L,L))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*( champ*(abs(y1-L/2)+abs(y2-L/2))/2 +inter*(y1-y2)**2) )
    return d
### Fonction qui calcule directement l'énergie libre
def distrib(ly,beta,h,j):
    tm = TransGauss(beta,h,j,ly)
    w,v = LA.eigh(tm)
    Z = np.sum(np.power(w,ly))
    ph = np.zeros(ly)
    for l in np.arange(ly):
        ph[l] = np.sum(np.power(w,ly)*np.power(v[l,:],2))
    return ph/Z,-1/np.log(w[-2]/w[-1])


###############################
# Mise en mémoire des données #
###############################
B = 0.01
T = 8

directory=sys.argv[1]
big_histo = np.loadtxt(directory)
print(big_histo)
fspace = [0]

sigma = np.empty((len(fspace),2))

fr = 10


h = big_histo[:,0]
LY = int(h[-1])+1
LY2 = int(LY/2)
ph = big_histo[:,1]/np.sum(big_histo[:,1])
plt.plot(h,ph,label="Monte Carlo")
popt,pcov= curve_fit(lambda h,m,l : np.log(p(h,alpha,m,l)), h[LY2-fr:LY2+fr],np.log(ph[LY2-fr:LY2+fr]),p0=[LY2,10])
plt.plot(h,p(h,alpha,*popt),'+',label='Fit Airy')
print(popt)
xiperp = abs(popt[1])
sigma[0] = T**2/(2*xiperp**3*B)
print(sigma)

Lmat = 200 ; lspace = np.arange(Lmat)
res = distrib(Lmat,1/T,B,1)
print(T,B)
ph = res[0]
plt.plot(lspace+LY2-Lmat/2,ph,label='Matrice de Transfert')

#plt.xlim([380,420])
#plt.ylim([1e-4,0.1])
#plt.yscale('log')
plt.xlabel('$h$')
plt.ylabel('$p(h)$')
#plt.plot(sigma[:,0],sigma[:,1])
plt.legend()
plt.savefig('airy-eq-kaw.pdf')
plt.show()
