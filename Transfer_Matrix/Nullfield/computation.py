#!/usr/bin/python
# -*- coding: utf8
import matplotlib.pyplot as plt
fsize=10
plt.rc("font",size=fsize)
plt.rc("font",family='serif')

import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import brentq
from numpy import linalg as LA
import sys
colors=['red','blue','green','rosybrown','sandybrown','gold','olivedrab']
nb_col = len(colors)
################################## 
### Déclaration fonctions ########
##################################
### Définition des matrices de transfert des différents modèles
def Trans(kbt,champ,inter,L):
    d=np.zeros((2*L+1,2*L+1))
    for y1,i in enumerate(d):
        for y2,j in enumerate(i):
            d[y1][y2] = np.exp(-kbt*inter*abs(y1-y2))
    return d
### Fonction qui calcule directement l'énergie libre
def intDiag(lx,ly,beta,h,j):
    tm = Trans(beta,h,j,ly)
    w,v = LA.eigh(tm) ; 
    return -np.log(max(w))/beta

def eigenvalue(beta,ly,prec):
    betexp = lambda beta : np.exp(-beta)
    theta = lambda eig,h,r : 1/np.tan(h*eig) - r*np.sin(eig)/(1-r*np.cos(eig))
    eigenvalue = lambda t,b : np.sinh(b)/(np.cosh(b)-np.cos(t))

    X = np.linspace(1e-20,1,pow(10,prec))
    func = theta(X[0],ly,betexp(BETA))
    for i,x in enumerate(X[1:]) :
        functmp = theta(x,ly,betexp(BETA))
        if(func*functmp < 0 ):
            zero = brentq(theta,X[i-1],x,args=(ly,betexp(beta)))
            return -np.log(eigenvalue(zero,beta))/beta
        else :
            func = functmp
    return False

############################
# Déclaration matrices
TC = 2/np.log(1.+np.power(2,0.5)) ; BETA = 1/TC;
LX = 70

Jn = 1
Hmax = 0
Lspace = np.arange(3,100)
Jspace = np.round(np.linspace(1,2,Jn),0)

####### Données #####
#### Calcul de la magnétisation via TM
def fitfreemag(x,surf,bulk,ex):
    return surf + 2*bulk*x + ex/x**2 

surface = np.empty(Jn)
bulk = np.empty(Jn)
exces = np.empty(Jn)

axfree = plt.subplot(111)
#axfit  = plt.subplot(212)

#######################
### Énergies libres ###
#######################
for nb,J in enumerate(Jspace):
    FreTM = np.empty([np.size(Lspace),2])
    FreDean = np.empty([np.size(Lspace),2])
    for numero,LY in enumerate(Lspace):
        FreTM[numero] = [LY,intDiag(LX,LY,BETA,Hmax,J)]
        FreDean[numero] = [LY,eigenvalue(BETA,LY+1,4)]

    axfree.plot(FreTM[:,0],FreTM[:,1],label='Diagonalisation')
    axfree.plot(FreDean[:,0],FreDean[:,1],label='Dean')
#    popw,pcov = curve_fit(fitfreemag,FreTM[:,0],FreTM[:,1],p0=[0.1,1,1])
#    axfree.plot(FreTM[:-1,0],fitfreemag(FreTM[:-1,0],*popw),'x',color=colors[nb%nb_col])

#    surface[nb] = popw[0]
#    bulk[nb] = popw[1]
#    exces[nb] = popw[2]


#axfree.set_yscale('log')
axfree.legend(loc='upper right')
axfree.set_xlabel('$L_Y$')
axfree.set_ylabel('$F(L_Y,J)$')

#######################
##### Coefficients ####
#######################
#def fitlin(x,a,b) :
#    return a*x+b
#def fitpow(x,a,p) :
#    return a*np.power(x,p)
#
#axfit.plot(Jspace,surface,label='Terme surface')
#popw,pcov = curve_fit(fitlin,Jspace,surface,p0=[0,1])
#print 'Surface :', popw
#
#axfit.plot(Jspace,bulk,label='Terme bulk')
#popw,pcov = curve_fit(fitlin,Jspace,bulk,p0=[0.1,0])
#print 'Bulk    :', popw
#
#axfit.plot(Jspace,exces,label='Terme exces')
#popw,pcov = curve_fit(fitpow,Jspace,exces,p0=[0.1,-1])
#axfit.plot(Jspace,fitpow(Jspace,*popw),'x',label='Fit exces')
#print 'Exces   :',popw
#print 'TC      :',TC
#print 'TC**2   :',TC**2
#print 'Beta    :',BETA
#print 'Beta**2 :',BETA**2
#
#axfit.set_xlabel('$J$')
#axfit.legend(loc='lower right')

plt.savefig('nullmag.pdf')
plt.show()
